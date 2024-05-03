# updated 220726: added option to change gap penalty for bowtie alignment
# Python 3.8
import sys
import os
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

args = sys.argv

r1 = args[1]
r2 = args[2]
reference = args[3]
sample_name = args[4]
gap_open_penalty = args[5] # default value for this with bowtie2 is 5. Suggest increasing for variants with lots of missense mutations

ref_name = reference.split(".")[0]
subfolder_name = "get_pileup_on_" + sample_name

# ======================================================================
# load modules on erisone
# modules = ["bowtie2", "trimgalore", "cutadapt", "fastqc", "samtools"]
# for mod in modules:
#    command = "module load " + mod
#    os.system(command)

# bowtie2 version 2.2.2
# trimgalore version 0.4.1
# cutadapt version 1.18
# fastqc version FastQC v0.10.1
# samtools version 1.4.1-1

# ====================================================================
# TRIM AND ALIGN
if not os.path.exists(subfolder_name + "/alignment.sam"):
    # make the bowtie index
    command = 'bowtie2-build ' + reference + " " + ref_name
    os.system(command)

    # make the subfolder
    command = "mkdir " + subfolder_name
    os.system(command)

    # trim the reads, put output into the subfolder
    command = "trim_galore --paired " + r1 + " " + r2 + " -o " + subfolder_name
    os.system(command)

    # map the reads
    r1_trimmed = subfolder_name + "/" + r1.split(".")[0] + "_val_1.fq.gz"
    r2_trimmed = subfolder_name + "/" + r2.split(".")[0] + "_val_2.fq.gz"
    command = 'bowtie2 -x ' + ref_name + ' -1 ' + r1_trimmed + ' -2 ' + r2_trimmed + " -S " + subfolder_name + "/alignment.sam --rg-id " + sample_name + " --rg SM:" + sample_name + " --rdg " + gap_open_penalty + ",3" + " --rfg " + gap_open_penalty + ",3" # note the --rg sample name is negessary or else later gatk will not know there is only one sample per sam file
    os.system(command)
    assert(os.path.exists(subfolder_name + "/alignment.sam"))
else:
     print("Skipping alignment. File already exists.")


# ==================================================================================
# MAKE COVERAGE PLOTS

# sort sam file and make coverage plot
command = "samtools sort -o " + subfolder_name + "/sorted_alignment.sam " + subfolder_name + "/alignment.sam"
os.system(command)
command = "samtools depth -a " + subfolder_name + "/sorted_alignment.sam > " + subfolder_name + "/coverage_table"
os.system(command)
# make plot of the coverage
coverage = pd.read_csv(subfolder_name + "/coverage_table", sep="	")
plt.figure(figsize=(15, 5))
sns.lineplot(x=coverage.iloc[:, 1], y=coverage.iloc[:, 2])
plt.xlabel("base position")
plt.ylabel("coverage (red=0, orange=25, green=50)")
plt.axhline(y=0, color='r', linestyle='-')
plt.axhline(y=25, color='orange', linestyle='-')
plt.axhline(y=50, color='green', linestyle='-')
plt.savefig(subfolder_name + "/coverage_plot")

# ======================================================================================
# get pileup file

if not os.path.exists(reference + ".fai"):
    command = "samtools faidx " + reference + " > " + reference + ".fai"
    os.system(command)
else:
    print(".fai file already exists")
if not os.path.exists(reference.split(".")[0] + ".dict"):
    command = "gatk-launch CreateSequenceDictionary -R " + reference
    os.system(command)
else:
    print(".dict file already exists")

command = "samtools mpileup " + subfolder_name + "/sorted_alignment.sam > " + subfolder_name + "/pileup -f " + reference
os.system(command)
