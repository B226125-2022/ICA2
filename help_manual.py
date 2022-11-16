import sys, os, subprocess
import re

#TEST for help manual specifically:
#esearch -db protein -query "pyruvate dehydrogenase[Protein Name] AND txid4890[organism] NOT partial[Properties]" | efetch -format fasta > fungi.fasta
#not partial[Property] removes partial entries

esearch -db protein -query "(pyruvate dehydrogenase[Protein Name] AND txid4890[organism] NOT partial[Properties]" | efetch -format fasta > fungi.fasta
#not partial[Property] #removes partial entries