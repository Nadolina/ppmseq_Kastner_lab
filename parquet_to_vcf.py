import pandas as pd
import pyarrow.parquet as pq
import numpy as np
import argparse
import os
import pysam.bcftools

## Parsing parquet name from command line input 
p=argparse.ArgumentParser()
p.add_argument("-s","--sample", help = "sample name", required=True)
p.add_argument("-u","--ug_sample", help = "Ultima Genomics sample name", required=True)
p.add_argument("-p","--parquet", help = "parquet file path", required=True)
args = p.parse_args()
if args.sample:
    print ("Sample to filter: % s" % args.sample)
    sample=str(args.sample)
if args.ug_sample:
    print ("Ultima Genomics sample name: % s" % args.ug_sample)
    ug_sample=str(args.ug_sample)
if args.parquet:
    print ("Parquet file path: % s" % args.parquet)
    parquet=str(args.parquet)

current_working_directory = os.getcwd()

p_path=parquet
header_path='/data/Kastner_PFS/scripts/pipelines/ppmseq/pq_to_vcf_header.txt' ## To be used as header for the output VCF file.
workdir=os.path.dirname(parquet)
outfile=os.path.basename(parquet).replace("parquet","vcf")

try:
    p=pq.ParquetFile(p_path)
    df=p.read().to_pandas()
except Exception as e:
    print(f"Error reading Parquet file: {e}")
    exit(1)

try:
    with open(header_path, 'r') as file:
        content = file.read()
except FileNotFoundError:
    print(f"Error: The file '{header_path}' was not found.")
    exit(1)
except Exception as e:
    print(f"An error occurred while reading the file: {e}")
    exit(1)

columns=["chromosome","position","identifier","reference","alternate","quality","filter",
         f"format_{ug_sample}_GT",f"format_{ug_sample}_DP",f"format_{ug_sample}_DP_FILT",f"format_{ug_sample}_VAF",
         f"format_{ug_sample}_RAW_VAF",f"format_{ug_sample}_AD",f"format_{ug_sample}_SNVQ",f"format_{ug_sample}_st",
         f"format_{ug_sample}_FILT"]

selectcols_df=df[columns]
renamed_df=selectcols_df.rename(columns={f"format_{ug_sample}_GT": "GT", f"format_{ug_sample}_DP": "DP", f"format_{ug_sample}_DP_FILT": "DP_FILT",
                  f"format_{ug_sample}_VAF": "VAF", f"format_{ug_sample}_RAW_VAF": "RAW_VAF", f"format_{ug_sample}_AD": "AD", f"format_{ug_sample}_SNVQ": "SNVQ", 
                  f"format_{ug_sample}_st": "ST",f"format_{ug_sample}_FILT":"FILT"})

with open(f'{workdir}/{outfile}', 'w') as vcf_file:
        vcf_file.write(content+"\n")
        vcf_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + sample + "\n")

        for row in renamed_df.itertuples():
            
            firststring="\t".join([str(row.chromosome), str(row.position), str(row.identifier[0]),str(row.reference), str(row.alternate[0]),
                                str(row.quality[0]), str(row.filter[0][0])])   

            info="."

            if row.VAF[0] < 0.8:
                GT="0/1"
            elif row.VAF[0] >= 0.8:
                GT="1/1"
            else:
                GT="NA"

            AD=str(row.AD[0][0]) + ',' + str(row.AD[0][1])

            NUM_MIXED=0
            for READ in row.ST:
                if READ == 'MIXED':
                    NUM_MIXED += 1
            NUM_SNVQ_GT60=0
            for SNVQ in row.SNVQ:
                if SNVQ >= 60:
                    NUM_SNVQ_GT60 += 1

            NUM_PASS=0
            NUM_SOFT_FILT=0
            for STATE in row.FILT:
                if STATE==1:
                    NUM_PASS+=1
                elif STATE==0:
                    NUM_SOFT_FILT+=1
            

            formatstring="GT:DP:DP_FILT:VAF:RAW_VAF:AD:NUM_MIXED:NUM_SNVQ_GT60:NUM_PASS:FRACTION_PASS"
            samplestring=":".join([GT, str(row.DP[0]), str(row.DP_FILT[0]), str(row.VAF[0]),
                                    str(row.RAW_VAF[0]),AD, str(NUM_MIXED), str(NUM_SNVQ_GT60),str(NUM_PASS),str(round(NUM_PASS/(NUM_PASS+NUM_SOFT_FILT),4))])

            vcf_file.write("\t".join([firststring, info, formatstring, samplestring]) + "\n")

vcf_file.close()

pysam.bcftools.norm("-O", "z", "-o", f"{workdir}/{outfile.replace('.vcf','.norm.vcf.gz')}", f"{workdir}/{outfile}", "-f", "/data/Kastner_PFS/references/HG38/Homo_sapiens_assembly38.fasta", catch_stdout=False)
