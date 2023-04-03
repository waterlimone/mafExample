##################################################################################################
# Author:  Cameron McIlvenna                                                                     #
# Date:    3/29/23                                                                               #
# Purpose: To read in and analyse pubicly available .maf file to determine statistical           #
#          significance between patients with KRAS gene and patients without it.                 #
##################################################################################################

from os import listdir, system
import subprocess
from shutil import rmtree
from time import sleep
import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
import time


# Merges all valid maf files into the first maf file
def merge(path, first):
    if(first == True):
        cmd = "tail -n +8 " + path + "/*.maf" + " > merged.maf"
        system(cmd)
    else:
        cmd = "tail -n +9 " + path + "/*.maf" + " >> merged.maf"
        system(cmd)


def remove_non_general():
    # Collect all the mafs directories for removing non general samples
    mafs = listdir("mafs")
    system("touch merged.maf")
    removed_files = 0
    is_first_maf = True  # Used to know if this is the first accepted maf so that the header is kept in the merged maf

    for maf in mafs:
        annotation_path = "mafs/" + maf + "/annotations.txt"
        maf_path = "mafs/" + maf
        annotations = pd.read_csv(annotation_path, sep="\t")
        categories = annotations.category.values.tolist()

        # If length of categories list is 1 otherwise remove maf directory
        if(len(categories) == 1):
            # If category is not "General" remove the maf directory
            if(categories[0] != "General"):
                removed_files += 1
                rmtree(maf_path)

            # If the category is "General" keep that directory
            else:
                # Unzip the maf file
                unzip_cmd = "gunzip " + maf_path + "/*.gz"
                system(unzip_cmd)

                # Invoke merge function
                merge(maf_path, is_first_maf)
                is_first_maf = False

        else:
            removed_files += 1
            rmtree(maf_path)

    print("Removed: " + str(removed_files) + " files")
    lines_after_merge = subprocess.check_output("cat merged.maf | wc -l", shell=True).decode("ascii")[0:-1]
    print("Lines after Merge: " + str(lines_after_merge))
    

# Deduplicating tumor barcodes for patients that have more than one
def deduplicate():
    merged = pd.read_csv("merged.maf", sep="\t", dtype="string")
    tumor_barcodes = merged.Tumor_Sample_Barcode.values.tolist()
    patients = []

    # Reads patient numbers from tumor barcodes
    for barcode in tumor_barcodes:
        patients.append(barcode[8:12])

    # set() used to rid duplicates and list() to convert set back to a list
    patients = list(set(patients))

    # For each patient find and remove duplicates with the exception of the last item
    for patient in patients:
        r = re.compile("[A-Z]{4}-.{2}-" + patient, re.IGNORECASE)

        # Filter tumor_barcodes by patient number using regex and remove duplicates with set(),
        # then convert back to list
        patient_dupes = list(set(list(
                filter(r.match, tumor_barcodes))))
        patient_dupes.sort() # Sort duplicates by alphabetical order

        # If the length of the dupes list is greater than 1 (there are duplicates),
        # then remove duplicates with the exception of the last one deduplicating the patient tumor_barcodes
        if(len(patient_dupes) > 1):
            patient_dupes = patient_dupes[0:-1]
            merged = merged[~merged.Tumor_Sample_Barcode.isin(patient_dupes)] # Merged where patient_dupes is not in merged

    print("Lines after Deduplication: " + str(merged.shape[0]))
    return merged

def dedup():
    merged = pd.read_csv("merged.maf", sep="\t", dtype="string")
    tumor_barcodes = merged.Tumor_Sample_Barcode.values.tolist()
    # Reads patient numbers from tumor barcodes
    patients = merged["Tumor_Sample_Barcode"].str[8:12].tolist()


    # set() used to rid duplicates and list() to convert set back to a list
    patients = list(set(patients))

    kept_codes = []
    # For each patient find and remove duplicates with the exception of the last item
    for patient in patients:
        r = re.compile("[A-Z]{4}-.{2}-" + patient, re.IGNORECASE)

        # Filter tumor_barcodes by patient number using regex and remove duplicates with set(),
        # then convert back to list
        patient_dupes = list(set(list(
                filter(r.match, tumor_barcodes))))
        patient_dupes.sort() # Sort duplicates by alphabetical order
        kept_codes.append(patient_dupes[-1])

    kept_codes = "|".join(kept_codes)
    merged = merged[merged["Tumor_Sample_Barcode"].str.match(kept_codes)]

    print("Lines after Deduplication: " + str(merged.shape[0]))
    return merged
    
        
def severe_mutation_enrichment(deduplicated):
    
    # Regex for PolyPhen matching "probably_damaging" or "possibly_damaging"
    rPoly = re.compile("probably_damaging|possibly_damaging")
    # Regex for SIFT matching "deleterious" or "deleterious_low_confidence"
    rSIFT = re.compile("deleterious|deleterious_low_confidence")

    # Filters PolyPhen and SIFT values using regex into lists
    poly_enriched = list(filter(rPoly.match, deduplicated.PolyPhen.dropna().values.tolist()))
    SIFT_enriched = list(filter(rSIFT.match, deduplicated.SIFT.dropna().values.tolist()))

    # Filters pandas dataframe from poly_enriched or SIFT_enriched
    deduplicated = deduplicated[(deduplicated.PolyPhen.isin(poly_enriched))|(deduplicated.SIFT.isin(SIFT_enriched))]
    print("Lines after Enrichment: " + str(deduplicated.shape[0]))
    return deduplicated

# Finds the most mutated genes counted only once per patient
def most_mutated_genes(enriched):
    enriched = pd.DataFrame(enriched)
    
    # Switched to groupby() for performance gains. 5 hundreths of a CPU second from 155!
    print(enriched.groupby("Hugo_Symbol", group_keys=False)["Tumor_Sample_Barcode"].nunique().sort_values(ascending=False)[0:5])

def analysis(enriched):
    # Read in CDR File
    cdr = pd.read_csv("TCGA-CDR-SupplementalTableS1.xlsx - TCGA-CDR.tsv", sep="\t")
    cdr["PFI_Bool"] = True # Make new column and set every value to true

    # If PFI value is 1.0 meaning cancer came back, set the corresponding bool to False to signify lower survival
    cdr.loc[cdr["PFI"] == 1.0, "PFI_Bool"] = False 
    
    # Gets patient barcodes as a list
    patient_codes = cdr.bcr_patient_barcode.values.tolist() 

    # Fixes issue with dataframe not having linting (It will still work but features like autocomplete will not)
    enriched = pd.DataFrame(enriched) 

    before_missing_patient_removal = enriched.shape[0] # Number of patients before removal of patients not in the cdr from the enriched set
    enriched = enriched[enriched.Tumor_Sample_Barcode.str.contains("|".join(patient_codes), regex=True)] # Filters out patients not in the cdr from enriched set
    patients_missing = before_missing_patient_removal - enriched.shape[0] # Get the number of patients missing from cdr file
    print("Patients Missing From CDR: " + str(patients_missing))

    # Filters enriched for only patients with KRAS then string slices the code to match up with cdr
    patients_with_KRAS = list(set(enriched.loc[enriched["Hugo_Symbol"] == "KRAS"]["Tumor_Sample_Barcode"].values.tolist()))
    patients_with_KRAS = [patients[0:12] for patients in patients_with_KRAS]

     # Filters enriched for only patients without KRAS then string slices the code to match up with cdr
    patients_withouth_KRAS = list(set(enriched.loc[enriched["Hugo_Symbol"] != "KRAS"]["Tumor_Sample_Barcode"].values.tolist()))
    patients_withouth_KRAS = [patients[0:12] for patients in patients_withouth_KRAS]

    # Filters the cdr using the patients with and without KRAS for analysis
    cdr_with_KRAS = cdr[cdr.bcr_patient_barcode.isin(patients_with_KRAS)]
    cdr_without_KRAS = cdr[cdr.bcr_patient_barcode.isin(patients_withouth_KRAS)]

    # Runs Kaplan Meier Estimator for with and without KRAS sets
    kmf_w = KaplanMeierFitter()
    kmf_wo = KaplanMeierFitter()
    kmf_w.fit(durations = cdr_with_KRAS["PFI.time"], event_observed = cdr_with_KRAS["PFI_Bool"], label = "With KRAS")
    kmf_wo.fit(durations = cdr_without_KRAS["PFI.time"], event_observed = cdr_without_KRAS["PFI_Bool"], label = "Without KRAS")

    # Runs logrank test
    logrank_test(cdr_with_KRAS["PFI.time"], cdr_without_KRAS["PFI.time"], 
                 event_observed_A=cdr_with_KRAS["PFI_Bool"], event_observed_B=cdr_without_KRAS["PFI_Bool"]).print_summary()

    # Plots the Kaplan Meier Curve
    kmf_w.plot()
    kmf_wo.plot()
    plt.xlabel("Days Since Cancer Diagnosis")
    plt.ylabel("Probability of Survival")
    plt.title("Survival Probabilities for Patients With and Without KRAS")
    plt.savefig("kaplanMeier", dpi=600)
    print("Saved Figure")
    

if(__name__ == "__main__"):
    st2 = time.time()
    st = time.process_time()
    remove_non_general()
    # deduplicated = deduplicate()    
    deduplicated = dedup()
    enriched = severe_mutation_enrichment(deduplicated)
    most_mutated_genes(enriched)
    analysis(enriched)
    print("CPU time:", time.process_time() - st)
    print("Actual time:", time.time() - st2)