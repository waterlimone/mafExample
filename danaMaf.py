from os import listdir, system
import subprocess
from shutil import rmtree
from time import sleep
import pandas as pd
import re


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
    is_first_maf = True  # Used to know if this is the first accepted maf so that the header is kept int he merged maf

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
    merged = pd.read_csv("merged.maf", sep="\t")
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

def most_mutated_genes(enriched):
    gene_dict = {}

    mutated_genes = list(set(enriched.Hugo_Symbol.values.tolist()))
    
    for gene in mutated_genes:
        patients_with_gene = enriched.loc[enriched["Hugo_Symbol"] == gene]["Tumor_Sample_Barcode"].values.tolist()
        patients_with_gene = list(set(
            [patients[8:12] for patients in patients_with_gene]))
        
        gene_dict[gene] = len(patients_with_gene)
    
    gene_dict = sorted(gene_dict.items(), key=lambda item: item[1], reverse=True)
    gene_dict = dict(gene_dict[0:5])
   
    print(gene_dict)

if(__name__ == "__main__"):
    remove_non_general()
    deduplicated = deduplicate()
    enriched = severe_mutation_enrichment(deduplicated)
    most_mutated_genes(enriched)

    

