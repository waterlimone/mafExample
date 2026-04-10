##################################################################################################
# Author:  Cameron McIlvenna                                                                     #
# Date:    3/12/26                                                                               #
# Purpose: To read in and analyse pubicly available .maf file to determine statistical           #
#          significance between patients with KRAS gene and patients without it.                 #
##################################################################################################
from time import time, process_time
import polars as pl

def enrich(parquet_file):
    st2 = time()
    st = process_time()
    print("Enrich Step")

    # Scan parquet file and enrich for PolyPhen scores of possibly_damaging | probably_damaging
    # or SIFT scores of deleterious | deleterious_low_confidence
    mafs = (
            pl.scan_parquet(parquet_file)
            .filter(
                pl.col("PolyPhen").str.contains_any(["possibly_damaging", "probably_damaging"])
                | pl.col("SIFT").str.contains_any(["deleterious", "deleterious_low_confidence"])
            )
            .sink_parquet("merged_deduplicated_enriched_mafs.parquet")
    )
    # Scans the metadata of the parquet file for the row count
    parquet_size = (
        pl.scan_parquet("merged_deduplicated_enriched_mafs.parquet")
        .select(pl.count()).collect().item()
    )
    
    print(f"\tRows Left:   {parquet_size}")
    print(f"\tCPU time:    {process_time() - st}")
    print(f"\tActual time: {time() - st2}")