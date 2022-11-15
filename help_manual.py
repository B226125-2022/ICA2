import sys, os, subprocess
import re

esearch -db protein -query "(pyruvate dehydrogenase[Protein Name] AND txid4890[organism] NOT partial[Properties]" | efetch -format fasta > fungi.fasta
#not partial[Property] #removes partial entries