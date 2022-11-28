#!usr/bin/env python3

import sys, os, subprocess
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
import matplotlib.image as mpimg

STRING_SPLIT_LENGTH = 70

#####################FUNCTIONS#######################################################################################################
def find_count(esearch_output): #Function to get Count
    count_regex = re.compile(r"<Count>(\d+)<\/Count>")
    for line in process_output.readlines():
        find_result = count_regex.findall(line)
        if len(find_result) > 0 and find_result[0].isdigit():
            return int(find_result[0])
    # if it doesn't find anything
    return 0

#get unique species_names
def unique_species_names(filename):
    result = []
    fasta_interesting_stuff_regex = re.compile(r"(>.+\[?.+\]?$)")
    f = open(filename, "r")
    for line in f:
        matches = fasta_interesting_stuff_regex.findall(line)
        if len(matches) > 0:
            result.append(matches[0].strip('[').strip(']'))
    return set(result)

def get_key_and_value(fasta_file):
    f = open(fasta_file, "r")
    file_as_string = f.read()
    # species_sequence = r'^(?P<strain_species>.+\]$)\n(?P<sequence>[A-Z\n]+)'
    species_sequence = r"(^>\[?.+\]?$)\n([A-Z\n]+)"
    match_group_tuples = re.findall(species_sequence, file_as_string, re.M)
    return dict([(k, v.replace('\n', '')) for k, v in match_group_tuples])

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

if esearch_count > 1000:
    print("More than 1000 sequences found, more specific search terms recommended.")
    sys.exit()

print("Working on it...")

"""
This has had a tendency to fail on DNS requests when the tool calls for CURL.
This may occur for the user as well, and is unfortunately beyond my control. 
This may be due to rate limiting, i.e. making a request too many times. 

curl: (6) Could not resolve host: eutils.ncbi.nlm.nih.gov
 ERROR:  curl command failed ( Fri 25 Nov 23:43:44 GMT 2022 ) with: 6
-X POST https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi [...]
 WARNING:  FAILURE ( Fri 25 Nov 23:43:44 GMT 2022 )
nquire -url https://eutils.ncbi.nlm.nih.gov/entrez/eutils/ efetch.fcgi [...]"""

cmd = "esearch -db protein -query \""+protein_fam_taxonomy+"\" | efetch -format fasta > "+fasta_file_name+""
process = subprocess.Popen(cmd, -1, shell=True, text=True, stdout=subprocess.PIPE)
process.wait()
process_output = process.stdout

print("Process complete, please check work directory if you would like to see the full fasta file output")

##############SEQUENCE PROCESSING#############################################################################################
my_file = open("sequence_general_info.txt","w")

while True:
    sequence_process = input("Would you like to continue with sequence processing? (y/n): ").lower()
    if sequence_process == "y" :
        species_names = unique_species_names(fasta_file_name)
        print("Unique species list of names is written to the file named 'sequence_general_info.txt'")        
        for names in species_names:
            my_file.write(names + "\n")
        break
    elif sequence_process == "n":
        print("Please feel free to try other protein families/taxonomy IDs")
        sys.exit()
    print("Please provide a y/n answer.\n")

my_file.close()

#fasta_file = input("Please input the fasta file name for further data processing: ")

species_sequence = get_key_and_value(fasta_file_name)


###############STATISTICS###################################################################################################
values = [len(value) for value in species_sequence.values()]
print("The average sequence length across species is: ", sum(values)/len(species_sequence))
print("The median sequence length is: ", np.median(values))

average = sum(values)/len(species_sequence)
median = np.median(values)

Q3, Q1 = np.percentile(values, [75 ,25])

print("The Q3 sequence length is: ", Q3)
print("The Q1 sequence length is: ", Q1)

IQR = Q3 - Q1
upper_boundary = Q3 + (1.5 * IQR)
print("Upper Boundary sequence length is: ", upper_boundary)
lower_boundary = Q1 - (1.5 * IQR)
print("Lower Boundary sequence length is: ", lower_boundary)

print("Outliers are defined as 1.5 multiplied by Q1 and Q3. If values are below and above these quartiles respectively, they will appear below.")

outliers = {}
for (k, v) in species_sequence.items():
    if len(v) < lower_boundary or len(v) > upper_boundary:
        print(f"Outlier {k} with length {len(v)}")
        outliers[k] = v

if len(outliers.keys()) == 0:
    print("No outliers have been found.")
else:
    delete_outlier = input("Would you like to delete these outliers? (y/n) (if you select n, the process will continue with the original fasta file): ").lower()
    if delete_outlier == "y":
        print(f"Attempting to remove {len(outliers.keys())} items")
        list(map(species_sequence.pop,outliers.keys()))
        print(f"Overwriting fasta file...")
        with open(fasta_file_name,"w") as write_fasta:
            for (key, value) in species_sequence.items():
                write_fasta.write(key)
                write_fasta.write("\n")
                values_list = [value[i:i+STRING_SPLIT_LENGTH] for i in range(0, len(value), STRING_SPLIT_LENGTH)]
                write_fasta.write("\n".join(values_list))
                write_fasta.write("\n")
            write_fasta.flush()
            write_fasta.close()
    else:
        print("No deletions, continuing to summary stats and alignment...")

##############DATAFRAME TO CSV###############################################################################################
df = pd.DataFrame( { "Protein Family" : protein_fam, 
"Taxonomy ID" : taxonomy, 
"Average Sequence Length" : average, 
"Quartile 1" : Q1, 
"Quartile 3" : Q3, 
"Interquartile Range" : IQR
}, index=[0])

csv_name = "summarised_stats_and_processed_data.csv"
df.to_csv(csv_name, sep=",", header=True)
print("Please check your directory to get the summarised data. File name: summarised_stats_and_processed_data.csv")

print("\nContinuing to sequence alignment, please wait a moment...")
##############ALIGNING SEQUENCES##############################################################################################
msf_file = f"{fasta_file_name[:-3]}.msf"

while True:
    if os.path.isfile(msf_file):
        print("You have specified a file name which exists. The msf file will be overwritten.")
        clustalo_cmd = f"clustalo -v -i {fasta_file_name} -o {fasta_file_name[:-3]}.msf --outfmt=msf --threads=20 --force"
        print("Aligning...")
        process = subprocess.Popen(clustalo_cmd, -1, shell=True, text=True, stdout=subprocess.PIPE)
        process.wait()
        process_output = process.stdout
        break
    else:
        clustalo_cmd = f"clustalo -v -i {fasta_file_name} -o {fasta_file_name[:-3]}.msf --outfmt=msf --threads=20"
        print("Aligning...")
        process = subprocess.Popen(clustalo_cmd, -1, shell=True, text=True, stdout=subprocess.PIPE)
        process.wait()
        process_output = process.stdout
    break

# clustalo_cmd = f"clustalo -v -i {fasta_file_name} -o {fasta_file_name[:-3]}.msf --outfmt=msf --threads=20"
# print("Aligning...")
# process = subprocess.Popen(clustalo_cmd, -1, shell=True, text=True, stdout=subprocess.PIPE)
# process.wait()
# process_output = process.stdout

output_file_infoalign = input("Please input a name for your sequence file for infoalignment: ")
infoalign_cmd = f"infoalign -sequence {msf_file} -outfile {output_file_infoalign}"

print("Processing with Infoalign... ")
process = subprocess.Popen(infoalign_cmd, -1, shell=True, text=True, stdout=subprocess.PIPE)
process.wait()
process_output = process.stdout

print(f"Infoalignment complete, please check the work directory for alignment information named : {output_file_infoalign} which contains general information of the sequence length, alignment length gaps, identity%, similarity% and more")

##############PLOTTING SEQUENCE CONSERVATION##################################################################################
while True:
    plotting = input("\nWould you like to continue to plotting? (y/n): ".lower())
    if plotting == "y":
        print("If you are ssh-ed into a server, make sure that you used the -Y option to make sure graphics work. Image loading may take a while...\nPlease close the graph to continue motif scanning")
        output_graph = f"{fasta_file_name[:-3]}"
        os.system(f"plotcon -sequence {msf_file} -winsize 4 -graph png -goutfile {output_graph}")
        img = mpimg.imread(f"{output_graph}.1.png")
        imgplot = plt.imshow(img)
        plt.axis("off")
        plt.show()
        break
    elif plotting == "n":
        print("Continuing to motif scanning")
        break
    print("Please provide a y/n answer.\n")

##############DETERMINING MOTIFS##############################################################################################
motif_output = "motifs"
total_outputs = []
value_filepath = "file_for_patmatmotifs.fa"

print("\nScanning for motifs...")

for (k, v) in species_sequence.items():
    with open(value_filepath, "w") as value_file:
        value_file.write(f"{k}\n{v}")
    subprocess.call(['bash', '-c', f"patmatmotifs -sequence {value_filepath} -outfile {motif_output}"])
    one_seq = subprocess.check_output(f'cat {motif_output}', shell=True).decode("UTF-8")
    
    with open("total_outputs.txt", "a") as total_outputs:
        total_outputs.write(one_seq)

print("Please check the file name 'total_outputs.txt' for motifs found across the fasta file.")

##############OTHER BIOLOGICAL INPUTS#########################################################################################
print("\nContinuing to other Biological Outputs of Interest...")
print("Getting peptide properties...")

one_species = "one_species.fa"

for index, (k, v) in enumerate(species_sequence.items()):
    print("Index: ", index, "\nSpecies: ", k)

select_species = input("Please select an index number to do further analysis on a specific species: ")
species_of_interest_key = list(species_sequence.keys())[int(select_species)]
species_of_interest_value = list(species_sequence.values())[int(select_species)]

with open(one_species, "w") as one_species_file:
    one_species_file.write(f"{species_of_interest_key}\n{species_of_interest_value}")
subprocess.call(['bash', '-c', f"pepinfo -sequence {one_species} -graph png -outfile {one_species}"])

print(f"Process complete, please check file names pepinfo.1.png, pepinfo.2.png, and {one_species} to get residue characteristics in both graphical and tabular formats")

while True:
    display_png = input("Would you like to display said .png files? (y/n): ".lower())
    if display_png == "y":
        print("If you are ssh-ed into a server, make sure that you used the -Y option to make sure graphics work. Image loading may take a while...\nPlease close the graph to continue to further analysis")
        img = mpimg.imread("pepinfo.1.png")
        imgplot = plt.imshow(img)
        plt.axis("off")
        plt.show()
        img2 = mpimg.imread("pepinfo.2.png")
        imgplot = plt.imshow(img2)
        plt.axis("off")
        plt.show()
        break
    elif display_png == "n":
        print("Program completed")
        sys.exit()
    print("Please provide a y/n answer.\n")
