#!/bin/sh 

#SBATCH --job-name=srsnv-df-filter
#SBATCH --cpus-per-task=25
#SBATCH --mem=250g
#SBATCH --time=24:00:00
#SBATCH --gres=lscratch:800
#SBATCH --output=logs/lowmem_%j.out
#SBATCH --error=logs/lowmem_%j.err

while getopts "s:c:h" option; do
   case $option in
	s) SAMPLE=$OPTARG ;;
	h) # display Help
         Help
         exit;;
   esac
done

set -euo pipefail
module load bcftools
module load singularity

## docker/singularity containers 
. /usr/local/current/singularity/app_conf/sing_binds
export SINGULARITY_BINDPATH="/data/$USER,/data/Kastner_PFS,/fdb,/lscratch/$SLURM_JOB_ID:/tmp"
SIF_SRSNV=/data/Kastner_PFS/ppmSeq/container/ugbio_srsnv_1.18.0.sif ## from https://hub.docker.com/r/ultimagenomics/featuremap
SIF_FEATUREMAP=/data/Kastner_PFS/ppmSeq/container/featuremap_master_9e07d7a.sif ## from https://hub.docker.com/r/ultimagenomics/featuremap

## sample and reference variables
BASE=${SAMPLE}
TRAINING_REGIONS=/data/Kastner_PFS/references/HG38/ultima_genomics/ug_rare_variant_hcr.Homo_sapiens_assembly38.interval_list.gz
TRAINING_REGIONS_INDEX=${TRAINING_REGIONS}.tbi
XGBOOST_PARAMS=/data/Kastner_PFS/references/HG38/ultima_genomics/250628.xgboost_model_params.json

# single_read_snv_params (ppmSeq template)
TP_TRAIN_SET_SIZE=1500000
FP_TRAIN_SET_SIZE=1500000
TP_OVERHEAD=10.0              # tp_train_set_size_sampling_overhead
MAX_VAF_FOR_FP=0.05
MIN_COV_FILTER=20             # min_coverage_filter
MAX_COV_FACTOR=2.0            # max_coverage_factor
RANDOM_SEED=0
NUM_FOLDS=3                   # num_CV_folds

MEAN_COVERAGE_FILE=${BASE}.mean_coverage.txt
MEAN_COVERAGE=$(cat "$MEAN_COVERAGE_FILE")
echo "Mean coverage: $MEAN_COVERAGE"
COVERAGE_CEIL=$(printf "%.0f" "$(echo "$MEAN_COVERAGE * $MAX_COV_FACTOR" | bc -l)")
echo "Coverage ceiling: $COVERAGE_CEIL"


# 3. Prepare RAW (negative / FP labeling set from raw featuremap)
#   a) Restrict to training regions
singularity exec ${SIF_SRSNV} \
    bcftools view ${SAMPLE}/${BASE}.raw.featuremap.vcf.gz -T ${TRAINING_REGIONS} -Oz -o ${SAMPLE}/${BASE}.raw.training_regions.vcf.gz
singularity exec ${SIF_SRSNV} \
    bcftools index -f -t ${SAMPLE}/${BASE}.raw.training_regions.vcf.gz

#   b) Convert to parquet
singularity exec ${SIF_SRSNV} \
    featuremap_to_dataframe \
    --input ${SAMPLE}/${BASE}.raw.training_regions.vcf.gz \
    --output ${SAMPLE}/${BASE}.raw.training_regions.parquet \
    --drop-format GT AD X_TCM 

#   c) Filter + label (RAW_VAF <= MAX_VAF_FOR_FP) + downsample to FP_TRAIN_SET_SIZE
singularity exec ${SIF_SRSNV} \
    filter_featuremap \
    --in  ${SAMPLE}/${BASE}.raw.training_regions.parquet \
    --out ${SAMPLE}/${BASE}.raw.filtered.parquet \
    --stats ${SAMPLE}/${BASE}.raw.stats.json \
    --filter name=coverage_ge_min:field=DP:op=ge:value=${MIN_COV_FILTER}:type=region \
    --filter name=coverage_le_max:field=DP:op=le:value=${COVERAGE_CEIL}:type=region \
    --filter name=mapq_ge_60:field=MAPQ:op=ge:value=60:type=quality \
    --filter name=no_adj_ref_diff:field=ADJ_REF_DIFF:op=eq:value=0:type=quality \
    --filter name=bcsq_gt_40:field=BCSQ:op=gt:value=40:type=quality \
    --filter name=edist_le_10:field=EDIST:op=lt:value=10:type=quality \
    --filter name=alt_hmer_lt_7:field=X_HMER_ALT:op=lt:value=7:type=quality \
    --filter name=low_vaf:field=RAW_VAF:op=le:value=${MAX_VAF_FOR_FP}:type=label \
    --downsample random:${FP_TRAIN_SET_SIZE}:${RANDOM_SEED}

# 4. Prepare RANDOM SAMPLE (positive / TP labeling set)
#   a) Restrict to training regions
singularity exec ${SIF_SRSNV} \
    bcftools view ${SAMPLE}/${BASE}.random_sample.featuremap.vcf.gz -T ${TRAINING_REGIONS} -Oz -o ${SAMPLE}/${BASE}.rs.training_regions.vcf.gz
singularity exec ${SIF_SRSNV} \
    bcftools index -t ${SAMPLE}/${BASE}.rs.training_regions.vcf.gz

#   b) Convert to parquet
singularity exec ${SIF_SRSNV} \
    featuremap_to_dataframe \
    --input ${SAMPLE}/${BASE}.rs.training_regions.vcf.gz \
    --output ${SAMPLE}/${BASE}.rs.training_regions.parquet \
    --drop-format GT AD X_TCM

#   c) Filter + label (REF == ALT) + downsample to TP_TRAIN_SET_SIZE
singularity exec ${SIF_SRSNV} \
    filter_featuremap \
    --in  ${SAMPLE}/${BASE}.rs.training_regions.parquet \
    --out ${SAMPLE}/${BASE}.rs.filtered.parquet \
    --stats ${SAMPLE}/${BASE}.rs.stats.json \
    --filter name=coverage_ge_min:field=DP:op=ge:value=${MIN_COV_FILTER}:type=region \
    --filter name=coverage_le_max:field=DP:op=le:value=${COVERAGE_CEIL}:type=region \
    --filter name=mapq_ge_60:field=MAPQ:op=ge:value=60:type=quality \
    --filter name=no_adj_ref_diff:field=ADJ_REF_DIFF:op=eq:value=0:type=quality \
    --filter name=bcsq_gt_40:field=BCSQ:op=gt:value=40:type=quality \
    --filter name=edist_le_10:field=EDIST:op=lt:value=10:type=quality \
    --filter name=alt_hmer_lt_7:field=X_HMER_ALT:op=lt:value=7:type=quality \
    --filter name=ref_eq_alt:field=REF:op=eq:value_field=ALT:type=label \
    --downsample random:${TP_TRAIN_SET_SIZE}:${RANDOM_SEED}

#   d) Filter + negative label (RAW_VAF <= MAX_VAF_FOR_FP) + downsample to TP_TRAIN_SET_SIZE
singularity exec ${SIF_SRSNV} \
    filter_featuremap \
    --in  ${SAMPLE}/${BASE}.rs.training_regions.parquet \
    --out ${SAMPLE}/${BASE}.rs_neg.filtered.parquet \
    --stats ${SAMPLE}/${BASE}.rs_neg.stats.json \
    --filter name=coverage_ge_min:field=DP:op=ge:value=${MIN_COV_FILTER}:type=region \
    --filter name=coverage_le_max:field=DP:op=le:value=${COVERAGE_CEIL}:type=region \
    --filter name=mapq_ge_60:field=MAPQ:op=ge:value=60:type=quality \
    --filter name=no_adj_ref_diff:field=ADJ_REF_DIFF:op=eq:value=0:type=quality \
    --filter name=bcsq_gt_40:field=BCSQ:op=gt:value=40:type=quality \
    --filter name=edist_le_10:field=EDIST:op=lt:value=10:type=quality \
    --filter name=alt_hmer_lt_7:field=X_HMER_ALT:op=lt:value=7:type=quality \
    --filter name=ref_ne_alt:field=REF:op=ne:value_field=ALT:type=label \
    --filter name=low_vaf:field=RAW_VAF:op=le:value=${MAX_VAF_FOR_FP}:type=label \
    --downsample random:${TP_TRAIN_SET_SIZE}:${RANDOM_SEED}
