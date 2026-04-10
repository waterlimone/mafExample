##################################################################################################
# Author:  Cameron McIlvenna                                                                     #
# Date:    3/12/26                                                                               #
# Purpose: To read in and analyse pubicly available .maf file to determine statistical           #
#          significance between patients with KRAS gene and patients without it.                 #
##################################################################################################
from time import time, process_time
from os import listdir
import polars as pl

def remove_non_general(maf_path):
    st2 = time()
    st = process_time()
    print("Merge Step")
    total = len(listdir(maf_path))
    
    # Scan the annotations and track the path to that folder.
    annotations = pl.scan_csv(  
                                f"{maf_path}/*/annotations.txt", 
                                separator="\t", 
                                has_header=True, 
                                include_file_paths="path"
    )

    # Filter and collect the paths of the annotation files that only have 1 category and is General
    filtered = (
                annotations
                .filter((pl.len().over("path") == 1) & (pl.col("category") == "General"))
                .with_columns(path = pl.col("path").str.replace("annotations.txt", "*maf.gz"))
                .select(["category","path"]).collect()
    )

    # Scan the maf files where the annotations were General only and stream the compressed file
    # to a parquet file.
    pl.scan_csv(
                    filtered["path"].to_list(),
                    separator="\t",
                    comment_prefix="#").sink_parquet("merged_mafs.parquet"
    )
    

    filtered_total = filtered.select(pl.len()).item()
    removed_files = total - filtered_total
    print(f"\tTotal Files:   {str(total)}")
    print(f"\tFiles Removed: {str(removed_files)}")
    print(f"\tFiles Left:    {str(filtered_total)}")
    print(f"\tCPU time:      {process_time() - st}")
    print(f"\tActual time:   {time() - st2}")