from os import listdir, system
import subprocess
from shutil import rmtree
from time import sleep
import pandas as pd
import re
import matplotlib.pyplot as plt
from sksurv.nonparametric import kaplan_meier_estimator


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
    print("Lines after enrichment: " + str(deduplicated.shape[0]))
    return deduplicated

# Finds the most mutated genes counted only once per patient
def most_mutated_genes(enriched):
    gene_dict = {}

    # Gets list of mutated genes and removes duplicates with set() and converts back to list
    mutated_genes = list(set(enriched.Hugo_Symbol.values.tolist())) 
    
    # For each gene pull the Tumor_Sample_Barcodes associated with that gene and slice for patient number,
    # then use set() to remove patient duplicates for that gene,
    # Finally add gene to gene_dict with the length of the patients_with_gene being the number of times this 
    # gene was seen. 
    for gene in mutated_genes:
        patients_with_gene = list(set(enriched.loc[enriched["Hugo_Symbol"] == gene]["Tumor_Sample_Barcode"].values.tolist()))
        patients_with_gene = [patients[8:12] for patients in patients_with_gene]
        
        gene_dict[gene] = len(patients_with_gene)
    
    # Sort the dictionary in descending order and take a slice of the first 5 values and convert back to a dict
    gene_dict = sorted(gene_dict.items(), key=lambda item: item[1], reverse=True)
    gene_dict = dict(gene_dict[0:5])
   
    print("Top 5 Most Mutated Genes: " + str(gene_dict))

def analysis(enriched):
    cdr = pd.read_csv("TCGA-CDR-SupplementalTableS1.xlsx - TCGA-CDR.tsv", sep="\t")
    cdr["PFI_Bool"] = True
    cdr.loc[cdr["PFI"] == 1.0, "PFI_Bool"] = False
    print(cdr["PFI"])
    print(cdr["PFI_Bool"])

    patient_codes = cdr.bcr_patient_barcode.values.tolist()
    enriched = pd.DataFrame(enriched)

    before_missing_patient_removal = enriched.shape[0]
    enriched = enriched[enriched.Tumor_Sample_Barcode.str.contains("|".join(patient_codes), regex=True)]
    patients_missing = before_missing_patient_removal - enriched.shape[0]
    print("Patients Missing From CDR: " + str(patients_missing))

    patients_with_KRAS = list(set(enriched.loc[enriched["Hugo_Symbol"] == "KRAS"]["Tumor_Sample_Barcode"].values.tolist()))
    patients_with_KRAS = [patients[0:12] for patients in patients_with_KRAS]
    patients_withouth_KRAS = list(set(enriched.loc[enriched["Hugo_Symbol"] != "KRAS"]["Tumor_Sample_Barcode"].values.tolist()))
    patients_withouth_KRAS = [patients[0:12] for patients in patients_withouth_KRAS]


    cdr_with_KRAS = cdr[cdr.bcr_patient_barcode.isin(patients_with_KRAS)]
    cdr_without_KRAS = cdr[cdr.bcr_patient_barcode.isin(patients_withouth_KRAS)]

    time_with, survival_prob_with = kaplan_meier_estimator(cdr_with_KRAS["PFI_Bool"], cdr_with_KRAS["PFI.time"])
    time_without, survival_prob_without = kaplan_meier_estimator(cdr_without_KRAS["PFI_Bool"], cdr_without_KRAS["PFI.time"])


    plt.step(time_with, survival_prob_with, where="post", label="With KRAS")
    plt.step(time_without, survival_prob_without, where="post", label="Without KRAS")
    plt.ylabel("Probability of Survival $\hat{S}(t)$")
    plt.xlabel("time $t$")
    plt.legend(loc="best")
    plt.show()
    

if(__name__ == "__main__"):
    # remove_non_general()
    deduplicated = deduplicate()    
    enriched = severe_mutation_enrichment(deduplicated)
    most_mutated_genes(enriched)
    analysis(enriched)