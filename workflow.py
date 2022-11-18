#!usr/bin/env python3

import sys, os, subprocess
import re #for count
import pandas as pd

#####################FUNCTIONS#######################################################################################################
def find_count(esearch_output): #Function to get Count
    count_regex = re.compile(r'<Count>(\d+)<\/Count>')
    for line in process_output.readlines():
        find_result = count_regex.findall(line)
        if len(find_result) > 0 and find_result[0].isdigit():
            return int(find_result[0])
    # if it doesn't find anything
    return 0

def unique_species_names(filename):
    result = []
    fasta_interesting_stuff_regex = re.compile(r">.+(\[.+\]$)")
    f = open(filename, 'r')
    for line in f:
        matches = fasta_interesting_stuff_regex.findall(line)
        if len(matches) > 0:
            result.append(matches[0].strip('[').strip(']'))
            # print(f"found this: {matches}")
    return set(result)



##############USER INPUT PROTEIN FAMILY AND TAXONOMY##################################################################################
protein_fam = input("Enter protein family: ")+"[Protein Name]"
taxonomy = input("Enter Taxonomy ID/Division\n(if taxonomy ID, type as txid##### ; Division example \"Ascomycota\"): ")+"[organism]"
print("Please note that partial proteins are not included")

protein_fam_taxonomy = protein_fam+" AND " +taxonomy+" NOT partial[Properties]"

fasta_file_spec_string = """
Enter the name of the output fasta file.
Please do not include spaces in the name.
Please end the filename with '.fa':
"""

while True:
    fasta_file_name = input(fasta_file_spec_string)

    if not fasta_file_name.endswith(".fa"):
        print("Please provide a filename which ends in '.fa'\n")
        continue

    if len(fasta_file_name) < 4:
        print("Please provide a filename before the '.fa' extension.\n")
        continue
    
    if os.path.isfile(fasta_file_name):
        print("You have specified a file name which exists. Would you like to delete this file?")
        while True:
            del_file = input("Would you like to delete the current file? (y/n):\n").lower()
            if del_file == "y":
                try:
                    os.remove(fasta_file_name)
                    break
                except Exception as e:
                    print("Encountered exception in deleting " + fasta_file_name)
                    print(e)
                    sys.exit(1)
            elif del_file == "n":
                print("Please specify a different filename.")
                break
            else:
                print("You did not provide a valid answer. Please try again.")
    break

cmd = "esearch -db protein -query \""+protein_fam_taxonomy+"\"" 

process = subprocess.Popen(cmd, -1, shell=True, text=True, stdout=subprocess.PIPE) #can't use os.system since that doesn't actually create a standard output
process.wait() 
process_output = process.stdout

esearch_count = find_count(process_output)

if not esearch_count > 0:
    print("No search results found, please check your inputs and try again.")
    sys.exit()

print("Working on it...")

cmd = "esearch -db protein -query \""+protein_fam_taxonomy+"\" | efetch -format fasta > "+fasta_file_name+""
process = subprocess.Popen(cmd, -1, shell=True, text=True, stdout=subprocess.PIPE)
process.wait()
process_output = process.stdout

print("process complete, please check work directory if you would like to see the full fasta file output")

##############SEQUENCE PROCESSING#############################################################################################
my_file = open("sequence_general_info.txt","w")

while True:
    sequence_process = input("Would you like to continue with sequence processing? (y/n): ").lower()
    if sequence_process == "y" :
        species_names = unique_species_names(fasta_file_name)
        print("Unique species list is written to the file named 'sequence_general_info'.txt")
        for names in species_names:
            my_file.write(names + "\n")
        break
    elif sequence_process == "n":
        print("Please feel free to try other protein families/taxonomy IDs")
        sys.exit()
    print("Please provide a y/n answer.\n")

my_file.close()
##############PLOTTING SEQUENCE CONSERVATION##################################################################################


##############DETERMINING MOTIFS##############################################################################################
#patmatmotifs	Scan a protein sequence with motifs from the PROSITE database

##############OTHER BIOLOGICAL INPUTS#########################################################################################