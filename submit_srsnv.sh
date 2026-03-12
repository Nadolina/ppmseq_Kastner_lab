#!/bin/sh

# submit script to order operations of Ultima Genomics SRSNV pipeline 
# SRSNV tools vary in their compute requirements, so I am breaking the pipeline into separate jobs to better manage resource allocation
Help() {
   echo "Usage: $0 -s <sample_name> -c <cram_file>"
   echo "Options:"
   echo "  -s    Sample name (required)"
   echo "  -c    CRAM file path (required)"
   echo "  -h    Display this help message"
}


while getopts "s:c:h" option; do
   case $option in
	s) SAMPLE=$OPTARG ;;
	c) CRAM=$OPTARG ;;
	h) # display Help
         Help
         exit;;
   esac
done

if [ -z "${SAMPLE:-}" ] || [ -z "${CRAM:-}" ]; then
    echo "Error: Missing required arguments"
    echo "Usage: $0 -s <sample_name> -c <cram_file>"
    exit 1
fi  

scriptpth="$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')"
sourcedir="$(dirname $scriptpth)"

## the snvfind job runs for a long time but needs minimal resources
snvfind_jid=$(sbatch ${sourcedir}/snvfind.sh -s "$SAMPLE" -c "$CRAM")

## the featuremap_to_dataframe step needs a lot of memory and cpu, dependent on snvfind completion
featuremap_jid=$(sbatch --dependency=afterok:${snvfind_jid} ${sourcedir}/featuremap_df_filter.sh -s "$SAMPLE")

## the training step (5) requires fewer resources so I am separating it from steps 3 and 4
sbatch --dependency=afterok:${featuremap_jid} ${sourcedir}/srsnv_train.sh -s "$SAMPLE"

