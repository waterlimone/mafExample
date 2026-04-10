##################################################################################################
# Author:  Cameron McIlvenna                                                                     #
# Date:    3/12/26                                                                               #
# Purpose: To read in and analyse pubicly available .maf file to determine statistical           #
#          significance between patients with KRAS gene and patients without it.                 #
##################################################################################################
from time import time, process_time
import polars as pl

def most_mutated_genes(parquet_file):
    st2 = time()
    st = process_time()
    print("Top 5 Mutated Genes")

    # Group by Hugo_Symbol and output the top 5 most mutated genes
    mafs = (
            pl.scan_parquet(parquet_file)
            .select("Hugo_Symbol")
            .group_by("Hugo_Symbol")
            .count()
            .sort("count", descending=True)
            .head(5)
            .collect()
    )

    print(mafs)
    print(f"\tCPU time:    {process_time() - st}")
    print(f"\tActual time: {time() - st2}")