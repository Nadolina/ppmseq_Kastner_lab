# SNV variant calling from ppmSeq data 

This README is divided into two parts: 1) Kastner lab version of SR-SNV from Ultimate Genomics, and 2) Filtering SR-SNV variants. The first part is all code from Ultima Genomics and uses Ultima Genomics docker image tools. I just needed to break the code into modules and make some minor additions to make it friendly to NIH's Biowulf environment. The second part includes all novel code. In this part we filter the SR-SNV VCF output for exome regions, MIXED start tags and high SNVQ scores. 

## 1. Kastner lab version of SR-SNV from Ultima Genomics
Modifying Ultima Genomics sr-snv pipeline to work on NIH's slurm HPC + companion scripts for filtering. Our workflow largely follows UG's recommendations in: https://github.com/Ultimagen/healthomics-workflows/blob/main/workflows/single_read_snv/howto-single-read-snv.md. 

### Requirements 

1. The singularity images ugbio_srsnv_1.16.1.sif and featuremap_master_883abc3.sif: https://hub.docker.com/r/ultimagenomics/
2. bcftools (the biowulf installation is fine, and you should not need to install anything)

### Inputs 

1. The sample ID. Note there is our assigned ID (usually 4 digits), and then the Ultima Genomics generated ID (which takes the form "BML[0-9][-0-9][0-9]").
2. The path to the cram file.
3. The directory containing your cram file must also containing a JSON file with the same prefix as your cram. This is a file we receive with the cram from NISC, and it contains stats necessary for this pipeline.

### Important outputs 

1. [SAMPLE].featuremap_df.parquet - likely positive SNV calls identified by the model 
2. [SAMPLE].featuremap.vcf.gz - the same content as the *.featuremap_df.parquet
3. [SAMPLE].report.html - contains a lot of useful diagnostic/QC info metrics

### Workflow 

The original pipeline is compiled into one script which can be found here: https://github.com/Ultimagen/healthomics-workflows/blob/main/workflows/single_read_snv/howto-single-read-snv.md. But I found variable computational requirements for the modules in this script. So, I broke the script into three modules according to computational allocations, and manage job submissions with submit_srsnv.sh. 

The __first step__, snvfind.sh, runs with low memory and CPU allocations. It generates coverage information, and produces a raw and sampled VCF that they call "featuremaps". 

The __second step__, featuremap_df_filter.sh:
* generates a parquet style dataframe from the training regions in the featuremap
* performs 'pre-filtering', which seems to be some basic QC
* assigns false positive, true positive and negative labels to the training set

In the __final script__, called srsnv_train.sh, these training sets are used to develop a model for distinguishing between true and artifact SNPs. This model is then applied back to the whole raw featuremap, resulting in a final featuremap parquet and an HTML report. 

Please refer to the Ultima Genomics github linked earlier for more information on these functions, as all functions are all re-purposed from their workflow. 

### Getting started 
```
sbatch /data/Kastner_PFS/scripts/pipelines/ppmseq/submit_srsnv.sh -s [SAMPLE ID] -c [PATH TO CRAM]
```

## 2. Filtering SR-SNV variants 

### Requirements 

1. mamba
2. use of a mamba environment with dependencies not available through biowulf modules (ppmseq_env.yml)

### Mamba on biowulf 

There are a couple softwares not available on biowulf that are required to run this filtering workflow. To workaround this, I have created a mamba environment that can be shared, but users will need their own mamba installation. Please refer to https://hpc.nih.gov/docs/diy_installation/conda.html for details, but to summarize:
1. start an sinteractive
2. load the mamba_install module on biowulf
3. run mamba_install to install conda
4. to use mamba, you need to source the install
```
[user@biowulf]$ sinteractive --mem=20g --gres=lscratch:20
[user@cn3444]$ module load mamba_install
[user@cn3444]$ mamba_install
...
[user@cn3444]$ source myconda
[user@cn3444]$ mamba --help
```

### The mamba "environment"

Environments are useful programming tools to keep software installations and versions separated for different tasks or workflows, and avoid installation conflicts. They also allow users to maintain consistency in collaborative projects. 

Now that we have mamba install, we can set up the environment required for the ppmseq filtering workflow. To start: 
1.  copy the YAML file to your working directory
2.  source your mamba your mamba install if you haven't already
3.  clone the environment from the YAML
4.  activate the env
```
source myconda
mamba env create -f ppmseq-env.yml
mamba activate ppmseqenv
```
**Troubleshooting note**:
If the "source myconda" does **not** work, you can try to run the following, which worked for a colleague. 
```
eval "$(mamba shell hook --shell bash)" 
```
And then, try "mamba activate ppmseqenv" again. 

### Getting started 
The only input required is the sample name. 
```
sbatch --mem=20g -c 8 --gres=lscratch:100 /data/Kastner_PFS/scripts/pipelines/ppmseq/run_SNVQ_filter.sh -s [SAMPLE ID]
sbatch --mem=20g -c 8 --gres=lscratch:100 /data/Kastner_PFS/scripts/pipelines/ppmseq/run_SNVQ_filter.sh -s 5447_2
```

### Summary of the workflow 

The goal of this workflow is to remove low quality reads from variant records in the outputs produced by SR-SNV. To summarize, we:
*  filter the *featuremap.vcf.gz to exome intervals
*  convert the exome VCF to a parquet file, for ease of use with pandas
*  expand the variant records to one-record-per-read
*  remove records/reads that do not meet one of the following filters:
    +  SNVQ > 55 & MIXED start-tag
    +  OR SNVQ >= 60 (regardless of start-tag type)
*  re-aggregate reads per-variant, and assign genotype based on VAF
    +  \> 0.8 homozygous (alternate)
    +  < 0.8 heterozygous
*  generate a VCF for the retained varaints of each of the two filters
*  normalize the variant records with pysam.bcftools 

### Outputs

The final filtered and normalized VCFs will have the naming scheme:
[SAMPLE].featuremap.exome.SNVQ60.norm.vcf.gz 
[SAMPLE].featuremap.exome.SNVQ55_stMIXED.norm.vcf.gz 


    
