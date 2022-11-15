#!usr/bin/env python3

import sys, os, subprocess
import re #for motifs ?

protein_fam = input("Enter protein family: ")+"[Protein Name]"
taxonomy = input("Enter Taxonomy ID/Division:\n(if taxonomy ID type as txid#### ; Division example Ascomycota): ")+"[organism]"
print("Please note that partial proteins are not included")

protein_fam_taxonomy = protein_fam+" AND " +taxonomy+" NOT partial[Properties]"
#print(protein_fam_taxonomy)

#TEST for help manual specifically:
#esearch -db protein -query "pyruvate dehydrogenase[Protein Name] AND txid4890[organism] NOT partial[Properties]" | efetch -format fasta > fungi.fasta
#not partial[Property] removes partial entries

cmd="esearch -db protein -query \""+protein_fam_taxonomy+"\"" | efetch -format fasta > test.fa" 


