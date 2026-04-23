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
1. generates a parquet style dataframe from the training regions in the featuremap
2. performs 'pre-filtering', which seems to be some basic QC
3. assigns false positive, true positive and negative labels to the training set

In the __final script__, called srsnv_train.sh, these training sets are used to develop a model for distinguishing between true and artifact SNPs. This model is then applied back to the whole raw featuremap, resulting in a final featuremap parquet and an HTML report. 

Please refer to the Ultima Genomics github linked earlier for more information on these functions, as all functions are all re-purposed from their workflow. 

### Getting started 
```
sbatch /data/Kastner_PFS/scripts/pipelines/ppmseq/submit_srsnv.sh -s [SAMPLE NAME] -c [PATH TO CRAM]
```

## 2. Filtering SR-SNV variants 

### Requirements 

1. use of the python environment /data/Kastner_PFS/scripts/pipelines/ppmSeq-venv
2. /data/Kastner_PFS/scripts/pipelines/ppmseq/pq_to_vcf_header.txt


### The python venv 

1. updated pyarrow because the biowulf pyarrow is version 6.0.1 but the latest (Feb 2026) is v. 23.0.1.
2. vcf2parquet (https://github.com/natir/vcf2parquet)




