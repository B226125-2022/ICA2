#!usr/bin/env python3

import sys, os, subprocess
import re #for motifs ?


def find_count(esearch_output):
    count_regex = re.compile(r'<Count>(\d+)<\/Count>')
    for line in process_output.readlines():
        find_result = count_regex.findall(line)
        if len(find_result) > 0 and find_result[0].isdigit():
            return int(find_result[0])
    # if it doesn't find anything
    return 0

def find_interesting_values(filename):
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
fasta_file_name = input("Enter the name of the output fasta file\nPlease do not include spaces in the name, addition of \".fa\" is unecessary: ")+".fa"

if os.path.isfile(fasta_file_name) == True: #if file name is not in current path:
    print("File", fasta_file_name, "already exists, please use a different name or delete the current fasta file in the working directory")
    del_file = input("Would you like to delete the current file (Y/N): ")
    if del_file == "Y" :
        os.remove(fasta_file_name)
        print("Please rerun the program to name the file again.")
    else:
        print("Please pick a different name.")
    sys.exit()
#print(protein_fam_taxonomy)

cmd = "esearch -db protein -query \""+protein_fam_taxonomy+"\"" 
# cmd2 = "esearch -db protein -query \""+protein_fam_taxonomy+"\" | efetch -format fasta > "+fasta_file_name+""

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


set_of_things = find_interesting_values(fasta_file_name)
for thing in set_of_things:
    print(thing + "\n")

# with open(fasta_file_name, 'r') as fasta_file:
#     fasta_file_contents = fasta_file.readlines()
#     fasta_interesting_stuff_regex = re.compile(r">.+(\[.+\]$)")
#     matches = fasta_interesting_stuff_regex.match(fasta_file_contents)
#     for match in matches:
#         print(match)


# pseudo code
#if intermediate file check passes, use original command:

# os.system(cmd2)
print("process complete, please check work directory for the output fasta file")

##############SEQUENCE PROCESSING#############################################################################################
#species_identification = re.search(>)

##############PLOTTING SEQUENCE CONSERVATION##################################################################################


##############DETERMINING MOTIFS##############################################################################################
#patmatmotifs	Scan a protein sequence with motifs from the PROSITE database

##############OTHER BIOLOGICAL INPUTS#########################################################################################