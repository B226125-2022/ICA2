#!usr/bin/env python3

import sys, os, subprocess
import re #for motifs ?

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
    sys.exit()
#print(protein_fam_taxonomy)

cmd="esearch -db protein -query \""+protein_fam_taxonomy+"\"" 
cmd2 = "esearch -db protein -query \""+protein_fam_taxonomy+"\" | efetch -format fasta > "+fasta_file_name+""

process = subprocess.Popen(cmd, -1, shell=True, text=True, stdout=subprocess.PIPE) #can't use os.system since that doesn't actually create a standard output
process.wait() 
process_output = process.stdout

with open('intermediate_file.txt', 'w') as the_file:
    the_file.writelines(process_output.readlines()) 
    the_file.flush()
    the_file.close()

#intermediate file check
# fh = open("intermediate_file.txt")
# check_data = fh.readlines()
# for line in enumerate(check_data, start = 1): #gives you value and index at the same time
#     print(type(line))
    # zero_count = re.search('<Count>0', line)
    # if zero_count:
    #     print("Your search retrieved zero results. Please check the protein family and/or taxonomy inputs again.")

#pseudo code
#if intermediate file check passes, use original command:
# os.sys(cmd2)
# print("Working on it...")
# print("process complete, please check work directory for the output fasta file")

##############PLOTTING SEQUENCE CONSERVATION##################################################################################

##############DETERMINING MOTIFS##############################################################################################

##############OTHER BIOLOGICAL INPUTS#########################################################################################