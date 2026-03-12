#!/bin/sh 

#SBATCH --job-name=srsnv-train
#SBATCH --cpus-per-task=4
#SBATCH --mem=12g
#SBATCH --time=12:00:00
#SBATCH --gres=lscratch:100
#SBATCH --output=logs/lowmem_%j.out
#SBATCH --error=logs/lowmem_%j.err

while getopts "s:h" option; do
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
RANDOM_SEED=0
NUM_FOLDS=3                   # num_CV_folds

MEAN_COVERAGE_FILE=${SAMPLE}/${BASE}.mean_coverage.txt
MEAN_COVERAGE=$(cat "$MEAN_COVERAGE_FILE")
echo "Mean coverage: $MEAN_COVERAGE"

# 5. Train (NUM_FOLDS=3 per ppmSeq template)
FEATURES="REF:ALT:X_PREV1:X_NEXT1:X_PREV2:X_NEXT2:X_PREV3:X_NEXT3:X_HMER_REF:X_HMER_ALT:BCSQ:BCSQCSS:RL:INDEX:REV:SCST:SCED:SMQ_BEFORE:SMQ_AFTER:tm:rq:st:et:EDIST:HAMDIST:HAMDIST_FILT:l1:l2:l3:l4:l5:l6:l7:q2:q3:q4:q5:q6"

singularity exec ${SIF_SRSNV} \
    srsnv_training \
    --positive ${SAMPLE}/${BASE}.rs.filtered.parquet \
    --negative ${SAMPLE}/${BASE}.raw.filtered.parquet \
    --stats-positive ${SAMPLE}/${BASE}.rs.stats.json \
    --stats-negative ${SAMPLE}/${BASE}.rs_neg.stats.json \
    --stats-featuremap ${SAMPLE}/${BASE}.raw.stats.json \
    --mean-coverage ${MEAN_COVERAGE} \
    --training-regions $TRAINING_REGIONS \
    --k-folds ${NUM_FOLDS} \
    --model-params $XGBOOST_PARAMS \
    --features $FEATURES \
    --basename $BASE \
    --output . \
    --random-seed ${RANDOM_SEED} \
    --verbose

# 6. Inference
mkdir -p ${SAMPLE}/model_files
cp ${SAMPLE}/${BASE}.model_fold_*.json ${SAMPLE}/model_files/
cp ${SAMPLE}/${BASE}.srsnv_metadata.json ${SAMPLE}/model_files/srsnv_metadata.json

singularity exec ${SIF_FEATUREMAP} \
    snvqual ${SAMPLE}/${BASE}.raw.featuremap.vcf.gz ${SAMPLE}/${BASE}.featuremap.vcf.gz model_files/srsnv_metadata.json -v
bcftools index -t ${SAMPLE}/${BASE}.featuremap.vcf.gz

# 7. Report
singularity exec ${SIF_SRSNV} \
    srsnv_report \
    --featuremap-df ${SAMPLE}/${BASE}.featuremap_df.parquet \
    --srsnv-metadata ${SAMPLE}/model_files/srsnv_metadata.json \
    --report-path . \
    --basename ${BASE} \
    --verbose
