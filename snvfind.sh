#!/bin/sh 

#SBATCH --job-name=srsnv-low
#SBATCH --cpus-per-task=4
#SBATCH --mem=6g
#SBATCH --time=24:00:00
#SBATCH --gres=lscratch:200
#SBATCH --output=logs/lowmem_%j.out
#SBATCH --error=logs/lowmem_%j.err

while getopts "s:c:h" option; do
   case $option in
	s) SAMPLE=$OPTARG ;;
	c) CRAM=$OPTARG ;;
	h) # display Help
         Help
         exit;;
   esac
done

set -euo pipefail
module load bcftools
module load singularity
ls -l /usr/local/current/singularity/app_conf/sing_binds
. /usr/local/current/singularity/app_conf/sing_binds ##following NIH recommendations on singularity use https://hpc.nih.gov/apps/singularity.html
export SINGULARITY_BINDPATH="/data/$USER,/data/Kastner_PFS,/fdb,/lscratch/$SLURM_JOB_ID:/tmp"
SIF_SRSNV=/data/Kastner_PFS/ppmSeq/container/ugbio_srsnv_1.18.0.sif ## from https://hub.docker.com/r/ultimagenomics/featuremap
SIF_FEATUREMAP=/data/Kastner_PFS/ppmSeq/container/featuremap_master_9e07d7a.sif ## from https://hub.docker.com/r/ultimagenomics/featuremap

BASE=${SAMPLE}
CRAM=${CRAM}
CRAM_INDEX=${CRAM}.crai
CRAM_PREFIX=$(basename "$CRAM" | sed 's/.cram//g')
SORTER_STATS=/data/Kastner_PFS/ppmSeq/2025/${BASE}/${CRAM_PREFIX}.json
REF=/data/Kastner_PFS/references/HG38/Homo_sapiens_assembly38.fasta
TRAINING_REGIONS=/data/Kastner_PFS/references/HG38/ultima_genomics/ug_rare_variant_hcr.Homo_sapiens_assembly38.interval_list.gz
TRAINING_REGIONS_INDEX=${TRAINING_REGIONS}.tbi
XGBOOST_PARAMS=/data/Kastner_PFS/references/HG38/ultima_genomics/250628.xgboost_model_params.json
BED=/data/Kastner_PFS/references/HG38/ultima_genomics/wgs_calling_regions.without_encode_blacklist.hg38.bed
TRINUC_FREQ=/data/Kastner_PFS/references/HG38/ultima_genomics/ref_trinuc_freq.csv  # Optional: reference genome trinucleotide frequency file for random sampling

if [ ! -f "${CRAM}" ] || [ ! -f "$CRAM_INDEX" ] || [ ! -f "$SORTER_STATS" ]; then
    echo "Error: One or more required files not found (CRAM, CRAM index, or sorter stats)"
    exit 1
fi

# single_read_snv_params (ppmSeq template)
TP_TRAIN_SET_SIZE=1500000
FP_TRAIN_SET_SIZE=1500000
TP_OVERHEAD=10.0              # tp_train_set_size_sampling_overhead
MAX_VAF_FOR_FP=0.05
MIN_COV_FILTER=20             # min_coverage_filter
MAX_COV_FACTOR=2.0            # max_coverage_factor
RANDOM_SEED=0
NUM_FOLDS=3                   # num_CV_folds


# Derived random sample size (ceil(tp_train_set_size * overhead))
RANDOM_SAMPLE_SIZE=$(( TP_TRAIN_SET_SIZE * 10 ))  # 15000000

mkdir -p ${SAMPLE}

# 1. Compute downsampling rate from sorter stats
TOTAL_ALIGNED_BASES=$(jq -re '.total_aligned_bases // .total_bases // error("missing total_aligned_bases")' "$SORTER_STATS")
DOWNSAMPLING_RATE=$(awk -v num=$RANDOM_SAMPLE_SIZE -v den=$TOTAL_ALIGNED_BASES 'BEGIN{printf "%.12f", num/den}')
echo "Downsampling rate: $DOWNSAMPLING_RATE"

MEAN_COVERAGE_FILE=${BASE}.mean_coverage.txt
singularity exec ${SIF_SRSNV} \
    sorter_stats_to_mean_coverage \
    --sorter-stats-json ${SORTER_STATS} \
    --output ${SAMPLE}/${MEAN_COVERAGE_FILE}

if [ ! -f "${SAMPLE}/${MEAN_COVERAGE_FILE}" ]; then
    echo "Error: Mean coverage file not found"
    exit 1
fi

MEAN_COVERAGE=$(cat "${SAMPLE}/${MEAN_COVERAGE_FILE}")
echo "Mean coverage: $MEAN_COVERAGE"
COVERAGE_CEIL=$(printf "%.0f" "$(echo "$MEAN_COVERAGE * $MAX_COV_FACTOR" | bc -l)")
echo "Coverage ceiling: $COVERAGE_CEIL"

# 2. snvfind (raw + random sample)
# FeatureMap params (ppmSeq): min_mapq=60 padding=5 score_limit=100 exclude_nan_scores=true include_dup_reads=true
# surrounding_quality_size=20 reference_context_size=3 keep_supplementary=false
# cram_tags_to_copy joined by commas below
CRAM_TAGS="tm:Z:A:AQ:AQZ:AZ:Q:QZ:Z,a3:i,rq:f,st:Z:MIXED:MINUS:PLUS:UNDETERMINED,et:Z:MIXED:MINUS:PLUS:UNDETERMINED,MI:Z,DS:i,sd:i,ed:i,l1:i,l2:i,l3:i,l4:i,l5:i,l6:i,l7:i,q2:i,q3:i,q4:i,q5:i,q6:i"

singularity exec ${SIF_FEATUREMAP} \
    snvfind ${CRAM} ${REF} \
    -o ${SAMPLE}/${BASE}.raw.featuremap.vcf.gz \
    -f ${SAMPLE}/${BASE}.random_sample.featuremap.vcf.gz,${DOWNSAMPLING_RATE}${TRINUC_FREQ} \
    -v \
    -p 5 -L 100 -n -d -Q 20 -r 3 -m 60 -c ${CRAM_TAGS} -b ${BED}

singularity exec ${SIF_SRSNV} \
    bcftools index -f -t ${SAMPLE}/${BASE}.raw.featuremap.vcf.gz
singularity exec ${SIF_SRSNV} \
    bcftools index -f -t ${SAMPLE}/${BASE}.random_sample.featuremap.vcf.gz


