#!/bin/sh 

#SBATCH --job-name=snvq-filter
#SBATCH --cpus-per-task=10
#SBATCH --mem=16g
#SBATCH --time=12:00:00
#SBATCH --gres=lscratch:100
#SBATCH --output=logs/snvq_filter_%j.out
#SBATCH --error=logs/snvq_filter_%j.err

module load bcftools

## for running the SNVQ filter script with the somatic conda environment

Help()
{
        # Display help 
        echo "To filter a VCF produced by the SRSNV pipeline, run this script with the following arguments:"
        echo "-s <sample_name> : Sample name to filter (required)"
        echo "-h : Display this help message"

}

while getopts "s:h" option; do
   case $option in
	s) SAMPLE=$OPTARG ;;
	h) # display Help
         Help
         exit;;
   esac
done

if [ -z "${SAMPLE:-}" ]; then
    echo "Error: Missing required arguments"
    echo "Usage: $0 -s <sample_name>"
    exit 1
fi  

OUTDIR=${PWD}/${SAMPLE}
if [ ! -d "$OUTDIR" ]; then
    echo "Creating output directory: $OUTDIR"
    mkdir -p "$OUTDIR"
fi

## Check for *featuremap.vcf.gz index 

if [ -f "${OUTDIR}/${SAMPLE}.featuremap.vcf.gz.csi" ]; then
    echo "Index found for ${OUTDIR}/${SAMPLE}.featuremap.vcf.gz - proceeding with extracting exome..."
else 
    echo "No index found for ${OUTDIR}/${SAMPLE}.featuremap.vcf.gz - generating index with bcftools..."
    bcftools index --threads ${SLURM_CPUS_PER_TASK} ${OUTDIR}/${SAMPLE}.featuremap.vcf.gz
fi


## Filter VCF to exomes -------------------

REGIONS=/data/Kastner_PFS/Sofia/References/010925_InMergedHg38Pad50_OutBL38curated_noALTdecoyContig.bed
bcftools view --threads ${SLURM_CPUS_PER_TASK} \
    -R ${REGIONS} \
    -Oz -o ${OUTDIR}/${SAMPLE}.featuremap.exomes.vcf.gz \
    --write-index ${OUTDIR}/${SAMPLE}.featuremap.vcf.gz

## Looking for ambiguous reference alleles in the VCF, which cause errors in the filter script. 
AMBIG_REF=$(bcftools query -f '%CHROM,%POS,%REF' ${OUTDIR}/${SAMPLE}.featuremap.exomes.vcf.gz | grep -e B -e K -e M -e N -e R -e S -e W -e Y -w) 
for AMBIG in ${AMBIG_REF}; do
    echo "Ambiguous reference allele found: $AMBIG, please confirm and remove manually, and then re-submit run_SNVQ_filter.sh."
done

## Converting the exome-limited VCF to a parquet file for filtering -------------------

source /data/Kastner_PFS/scripts/pipelines/ppmSeq-venv ## required updated installations of vcf2paruqte and pyarrow for the parquet conversion to work.

vcf2parquet --input ${OUTDIR}/${SAMPLE}.featuremap.exomes.vcf.gz convert --output ${OUTDIR}/${SAMPLE}.featuremap.exomes.parquet


## Run filter -------------------

## ultima genomics sample name is different from our assigned sample name, so we need to extract it from the VCF file and provide it to the filter script
UG_SAMPLE=$(bcftools query -l ${OUTDIR}/${SAMPLE}.featuremap.exomes.vcf.gz | head -n 1)
echo "Ultima Genomics sample name extracted from VCF: ${UG_SAMPLE}"

python /data/Kastner_PFS/scripts/pipelines/ppmseq/SNVQ_filter.py -s ${SAMPLE} -u "${UG_SAMPLE}"