#!/usr/bin/env python3

import pandas as pd
import pyarrow.parquet as pq
import numpy as np
import argparse
import pysam.bcftools

## Parsing sample ID from command line input 
p=argparse.ArgumentParser()
p.add_argument("-s","--sample", help = "Pass sample name to filter script.", required=True)
p.add_argument("-u","--ug_sample", help = "Pass the sample name used by ultima genomics to filter script.", required=True)
args = p.parse_args()
if args.sample:
    print ("Sample to filter: % s" % args.sample)
    sample=str(args.sample)
if args.ug_sample:
    print ("Ultima Genomics sample name: % s" % args.ug_sample)
    ug_sample=str(args.ug_sample)

## Getting working directory
workdir=os.getcwd()

p_path=f'{workdir}/{sample}/{sample}.featuremap.exomes.parquet'
p=pq.ParquetFile(p_path)

## VARIABLES --------------

batch_size=100000
columns_to_remove=['SCST','SCED','MAPQ','EDIST','HAMDIST','HAMDIST_FILT', 'FILT_BITMAP','MQUAL','tm','a3','MI','DS','sd','ed','l1','l2',
                    'l3','l4','l5','l6','l7','q2','q3','q4','q5','q6','SMQ_BEFORE','SMQ_AFTER','ADJ_REF_DIFF']
columns_to_explode=['BCSQ','BCSQCSS','RL','INDEX','RN','DUP','REV','SNVQ','st','et','rq','FILT']

columns_to_remove_full=[f'format_{ug_sample}_{col}' for col in columns_to_remove]
columns_to_explode_full=[f'format_{ug_sample}_{col}' for col in columns_to_explode]
    
idxcols=["chromosome","position","reference","alternate"]
output_parquet1=f'{workdir}/{sample}/{sample}.featuremap.exomes.SNVQ55_stMIXED.parquet'
output_parquet2=f'{workdir}/{sample}/{sample}.featuremap.exomes.SNVQ60.parquet'

## FUNCTIONS ------------------------------------------------------------------

def explode_filter_batch(df):
    df.drop(columns=columns_to_remove_full, inplace=True)
    exploded_df=df.explode(columns_to_explode_full)
    filtered_df1=exploded_df[(exploded_df[f'format_{ug_sample}_SNVQ']>55) & (exploded_df[f'format_{ug_sample}_st']=='MIXED')]
    filtered_df2=exploded_df[exploded_df[f'format_{ug_sample}_SNVQ']>=60]
    return filtered_df1, filtered_df2

def reaggregate_batch(df):
    collist=[item for item in df.columns.tolist() if item not in idxcols]

    df=df.copy()
    df["record"] = list(
        zip(*(df[col] for col in collist))
    )

    ## the alternate column serves an index, along with chromosome, position and reference, for grouping to ensure all records are grouped correctly
    ## but alternate is an array, which is unhashable, so we need to convert it to a tuple for grouping
    df["alternate"] = df["alternate"].map(
        lambda x: tuple(x) if isinstance(x, np.ndarray) else x
    )       
    zipped_df = (df.groupby(idxcols, sort=False, as_index=False).agg(record=("record", list)))
    return zipped_df, collist

## To ensure the corresponding values are grouped in the correct order after re-aggregation, we created a "records" column that zips all the relevant columns into a tuple. 
## I want to now unzip these records back into their original columns.
def unzip_records(recs, n):
    # recs: list[tuple] (length may be 0)
    if not recs:
        return tuple([[] for _ in range(n)])
    return tuple(map(list, zip(*recs)))  # -> (list_a, list_b, list_c, ...)


## MAIN ---------------------------------------------------------------------

## initializing list to collect filtered batches, to be concatenated after filtering
## 2 lists because 2 filtering approaches, resulting in two final parquets
processed_batches1=[]
processed_batches2=[]
for i, batch in enumerate(p.iter_batches(batch_size=batch_size)):
    print (f"Processing batch {i+1}...")

    ## Filtering the batch according to two filters, 1) SNVQ > 55 & start-tag is MIXED (batch_filtered1), or 2) SNVQ >= 60 (batch_filtered2)
    batch_df=batch.to_pandas()
    batch_filtered1, batch_filtered2 = explode_filter_batch(batch_df)

    ## re-aggregating variant records in batch with SNVQ & MIXED start-tag filters ---

    batch_zipped1, collist1=reaggregate_batch(batch_filtered1)

    n = len(collist1)

    batch_zipped1[collist1] = pd.DataFrame(
        batch_zipped1["record"].apply(lambda recs: unzip_records(recs, n)).tolist(),
        columns=collist1,
        index=batch_zipped1.index
    )
    batch_zipped1.drop(columns=["record"], inplace=True)

    processed_batches1.append(batch_zipped1)

    ## re-aggregating variant records in batch with just SNVQ >= 60 filter -----

    batch_zipped2, collist2=reaggregate_batch(batch_filtered2)

    batch_zipped2[collist2] = pd.DataFrame(
        batch_zipped2["record"].apply(lambda recs: unzip_records(recs, n)).to_list(),
        columns=collist2,
        index=batch_zipped2.index
    )
    batch_zipped2.drop(columns=["record"],inplace=True)
    processed_batches2.append(batch_zipped2)

    del batch_df 

processed_df1=pd.concat(processed_batches1, ignore_index=True)
processed_df2=pd.concat(processed_batches2, ignore_index=True)

row_count = processed_df1.shape[0]
print(f"Number of rows after filtering with SNVQ > 55 AND start-tag=MIXED: {row_count}") 
row_count = processed_df2.shape[0]
print(f"Number of rows after filtering with SNVQ >= 60: {row_count}") 

processed_df1.to_parquet(output_parquet1, engine='pyarrow')
processed_df2.to_parquet(output_parquet2, engine='pyarrow')


