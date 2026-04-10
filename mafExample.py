##################################################################################################
# Author:  Cameron McIlvenna                                                                     #
# Date:    3/12/26                                                                               #
# Purpose: To read in and analyse pubicly available .maf file to determine statistical           #
#          significance between patients with KRAS gene and patients without it.                 #
##################################################################################################

from remove_non_general import remove_non_general
from deduplicate import deduplicate
from enrich import enrich
from most_mutated_genes import most_mutated_genes
from time import time, process_time

if(__name__ == "__main__"):
    st2 = time()
    st = process_time()

    remove_non_general("mafs")
    deduplicate("merged_mafs.parquet")
    enrich("merged_deduplicated_mafs.parquet")
    most_mutated_genes("merged_deduplicated_enriched_mafs.parquet")

    print(f"CPU time:    {process_time() - st}")
    print(f"Actual time: {time() - st2}")