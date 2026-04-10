##################################################################################################
# Author:  Cameron McIlvenna                                                                     #
# Date:    3/12/26                                                                               #
# Purpose: To read in and analyse pubicly available .maf file to determine statistical           #
#          significance between patients with KRAS gene and patients without it.                 #
##################################################################################################

from time import time, process_time
import polars as pl

def deduplicate(parquet_file):
    st2 = time()
    st = process_time()
    print("Deduplication Step")
    
    # Scans the metadata of the parquet file for the row count
    parquet_size = pl.scan_parquet(parquet_file).select(pl.count()).collect().item()
    print(f"\tRows Before: {parquet_size}")

    # Scan the parquet file that was merged, create a column of the patient ids to have unifying group.
    # Using a stable sort make sure the Tumor_Samples are still in alphabetical order, and then sort by
    # unifying patient group as the primary sort. Finally filter for the max Tumor_Samples over the window
    # of unifying patient codes grabbing the whole group. Uses sink_parquet to keep memory pressure down
    # and for NextFlow integration writing this step to a file.
    mafs = (
            pl.scan_parquet(parquet_file)
            .with_columns((
                pl.col("Tumor_Sample_Barcode")
                .str.slice(8,4)
                .alias("Patient_Codes"))
            )
            .sort("Tumor_Sample_Barcode")
            .sort("Patient_Codes")
            .filter(
                pl.col("Tumor_Sample_Barcode") == pl.col("Tumor_Sample_Barcode")
                .max()
                .over("Patient_Codes")
            )
            .sink_parquet("merged_deduplicated_mafs.parquet")
    )

    # Scans the metadata of the parquet file for the row count
    parquet_size = pl.scan_parquet("merged_deduplicated_mafs.parquet").select(pl.count()).collect().item()

    print(f"\tRows Left:   {parquet_size}")
    print(f"\tCPU time:    {process_time() - st}")
    print(f"\tActual time: {time() - st2}")