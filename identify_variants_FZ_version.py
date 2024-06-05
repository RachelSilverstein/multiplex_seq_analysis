# go into each folder whose name starts with 'call_variants_on' in the 'main_dir' directory and see what variant is in the pileup file
# summarize the pileup file in a 'pileup_summary.csv'
# output a file on the main directory called 'output_variants.csv' which contains a list of all the variant identities and their folder names

import pandas as pd
import os
from Bio.Seq import Seq

main_dir = "."
MAX_INSERTION_FREQ = 0.2  # call a variant an insertion when frequency of reads with deletions above this fraction
MAX_DELETION_FREQ = 0.2  # call a variant a deletion when frequency of reads with deletions above this fraction
call_variant_when_identity_below = 0.95  # include a position in the pileup summary table if the fraction of bases which match the reference at that position is below this fraction
include_alleles_above_freq = 0.01 # include all the alleles at a position in the pileup summary table that occur above this frequency
MIXED_VARIANT_FRACTION = 0.05 # maximim fraction of impurities at one of the mutated positions without calling the variant mixed

ALLOWED_PILEUP_CHARACTERS = ["A", "C", "G", "T", ".", "+", "*", "-", "1", "2", "3", "4", "5", "6", "7", "8", "9", "0"]
mutated_positions = [3441, 3442, 3443, 3444, 3445, 3446, 3690, 3691, 3692, 3693, 3694, 3695, 4035, 4036, 4037, 4041, 4042, 4043, 4047, 4048, 4049] # NOTE THESE BASES ARE SPECIFICALLY FOR REFERENCE FZ-Cas9_RTW3027_extraction
recoded_positions = {}  # example {4106: "G", 4112: "T"}


## ---------------------------------------------------------------------------------------------------------------------------------
dir_contents = os.listdir(main_dir)

subdirs = []
for element in dir_contents:
    if "call_variants_on_" in element:
        subdirs.append(element)

samples = []
var_names = []
D1135s = []
S1136s = []
G1218s = []
E1219s = []
R1333s = []
R1335s = []
T1337s = []
low_reads_warnings = []
min_read_positions = []
var_types = []
coding_variant_warnings = []
coding_variant_over_50 = []

def calculate_average_q_score(string):
    """Take a string of q score characters and calculate the average"""
    score_sum = 0
    for char in string:
        score_sum += ord(char) - 33
    return score_sum / len(string)

def pileup_to_list(bases):
    """Takes a row of bases from mpileup file and converts it to a list of alleles.
    Helper for get_mismatch table"""
    bases = bases.replace(",", ".")
    bases = bases.replace("a", "A")
    bases = bases.replace("c", "C")
    bases = bases.replace("g", "G")
    bases = bases.replace("t", "T")
    bases_list = []
    i = 0
    while i < len(bases):
        char = bases[i]
        if char in [".", "A", "T", "C", "G", "*"]:
            bases_list.append(char)
            i += 1
        elif char == "^":  # start of a read, remove the quality score which follows as this can be confused with indel length
            i += 2
        elif char == "+" or char == "-":
            num_digits = 0
            j = 1
            while j < 5 and i + j < len(bases):
                if bases[i + j] in ["1", "2", "3", "4", "5", "6", "7", "8", "9", "0"]:
                    num_digits = j
                    j += 1
                else:
                    break
            indel_length = int(bases[i+1:i+num_digits+1])
            indel_string = bases[i:i+indel_length + num_digits + 1]
            i += indel_length + num_digits + 1
            bases_list.append(indel_string)
        elif char in ["$"]:
            i += 1
            continue
        else:
            print("unrecognized character")
            print(char)
            raise Exception()
    return bases_list


def get_mismatch_table(pileup, call_variant_when_identity_below, include_alleles_above_freq):
    """Get a list of mismatches in a variant from the pileup file.
    Create a row in the mismatch table for each position that has less than 90% identity with the reference or
    an indel following that position at greater than 10% of coverage.
    For each position listed, list all alleles occurring in over 5% of reads"""
    positions = []
    ref_nt = []
    alt_nt = []
    alt_frequencies = []
    av_q_scores = []
    for i in range(len(pileup)):
        row = pileup.iloc[i, ]
        bases = row["bases"]
        bases = pileup_to_list(bases)
        matches = bases.count(".")
        av_q_score = calculate_average_q_score(row["quality"])
        if row["coverage"] != 0:
            match_freq = matches/row["coverage"]
            next_insertion_freq = sum([element[0] == "+" for element in bases])/row["coverage"]
            next_deletion_freq = sum([element[0] == "-" for element in bases])/row["coverage"]
            if match_freq < call_variant_when_identity_below or \
                    next_insertion_freq > (1-call_variant_when_identity_below) or \
                    next_deletion_freq > (1-call_variant_when_identity_below):
                alt_alleles = []
                frequencies = []
                for allele in set(bases) - {"."}:
                    frequency = bases.count(allele)/row["coverage"]
                    if frequency > include_alleles_above_freq:
                        frequencies.append(frequency)
                        alt_alleles.append(allele)
                positions.append(row["position"])
                ref_nt.append(row["ref_nt"])
                alt_nt.append(alt_alleles)
                alt_frequencies.append(frequencies)
                av_q_scores.append(av_q_score)
    output = pd.DataFrame({"position": positions, "ref_nt": ref_nt, "alt_nt": alt_nt, "alt_frequencies": alt_frequencies, "av_q_score": av_q_scores})
    return output


def get_variant_type(mismatch_table,
                     max_insertion_freq,
                     max_deletion_freq,
                     mixed_variant_fraction):
    """Get the type of variant based on the mismatch table.
    Can be either 'missense or nonsense', 'insertion', 'deletion' 'insertion and deletion' or 'mixed variant'"""
    var_type = "missense or nonsense"
    insertion = False
    deletion = False
    mixed_variant = False
    for i in range(len(mismatch_table)):
        row = mismatch_table.iloc[i, ]
        position = row["position"]
        alt_nt = row["alt_nt"]
        alt_freqs = row["alt_frequencies"]
        if any([allele[0] == "+" for allele in alt_nt]):
            indexes = [i for i in range(len(alt_nt)) if alt_nt[i][0] == "+"]
            freqs = [alt_freqs[index] for index in indexes]
            if any([freq > max_insertion_freq for freq in freqs]):
                insertion = True
        if any([allele[0] == "-" for allele in alt_nt]):
            indexes = [i for i in range(len(alt_nt)) if alt_nt[i][0] == "-"]
            freqs = [alt_freqs[index] for index in indexes]
            if any([freq > max_deletion_freq for freq in freqs]):
                deletion = True
        if "*" in alt_nt:
            indexes = [i for i in range(len(alt_nt)) if alt_nt[i] == "*"]
            freqs = [alt_freqs[index] for index in indexes]
            if any([freq > max_deletion_freq for freq in freqs]):
                deletion = True
        if (position in mutated_positions) and (insertion is False) and (deletion is False):
            pure_alt_nt = any([freq > (1-mixed_variant_fraction) for freq in alt_freqs])  # there is a single alternate nucleotide taking up most of reads
            pure_wt_nt = sum([freq for freq in alt_freqs]) < mixed_variant_fraction  # the sum of the impurities frequencies is below the allowed threshhold
            if not (pure_alt_nt or pure_wt_nt):
                mixed_variant = True
    if insertion and deletion:
        var_type = "deletion and insertion"
    elif insertion:
        var_type = "insertion"
    elif deletion:
        var_type = "deletion"
    elif mixed_variant:
        var_type = "mixed variant"
    return var_type

def get_new_sample_name(old_name):
    """ if the samples are not named with the zero in names such as A01 A02 (instead A1 A2 etc.)
    Then add the zeros so that excel will sort the names correctly"""
    for letter in ["A", "B", "C", "D", "E", "F", "G", "H"]:
        for num in ["1", "2", "3", "4", "5", "6", "7", "8", "9"]:
            well = letter + num + "_"
            if well in old_name:
                return old_name.replace(well, letter + "0" + num + "_")
    return old_name

for dir in subdirs:
    # get the name of the sample from the directory name
    sample = dir.split("_")[3:]
    sample = "_".join(sample)
    new_sample_name = get_new_sample_name(sample)
    print(new_sample_name)
    samples.append(new_sample_name)
    # make sure the read counts are ok
    coverage = pd.read_csv(main_dir + "/" + dir + "/coverage_table", sep="	")
    coverage = coverage.iloc[17:5222, 2]
    min_coverage = min(coverage)
    if min(coverage) < 25:
        low_reads_warning = "True - min coverage is " + str(min_coverage)
    else:
        low_reads_warning = "False"
    # output the position of the lowest coverage
    minimum_index = coverage.idxmin(axis=0)
    min_read_positions.append(minimum_index)
    low_reads_warnings.append(low_reads_warning)

    # get pileup table
    pileup = pd.read_csv(main_dir + "/" + dir + "/pileup", sep="	",
                         names=["ref_name", "position", "ref_nt", "coverage", "bases", "quality"])
    # make mismatch table
    mismatch_table = get_mismatch_table(pileup, call_variant_when_identity_below, include_alleles_above_freq)
    mismatch_table.to_csv(main_dir + "/" + dir + "/pileup_summary.csv")





    # get the summary data on the variant and add to the table
    var_type = get_variant_type(mismatch_table, MAX_INSERTION_FREQ, MAX_DELETION_FREQ, MIXED_VARIANT_FRACTION)
    var_types.append(var_type)
    print(var_type)


    if var_type == "missense or nonsense" and min_coverage > 0:
        nts_at_allowed_positions = []
        # originally set all the positions to be their WT nucleotide:
        for position in mutated_positions:
            nt = pileup[pileup["position"] == position]["ref_nt"].values
            nts_at_allowed_positions.extend(nt)

        # mutate the positions if they are mutated in the pileup summary file
        for i in range(len(mismatch_table)):
            row = mismatch_table.iloc[i, ]
            position = row["position"]
            alt_alleles = row["alt_nt"]
            alt_freqs = row["alt_frequencies"]
            if position in mutated_positions:
                max_index = alt_freqs.index(max(alt_freqs))
                new_allele = alt_alleles[max_index]
                assert(len(new_allele) == 1) # should only be subbing in 1 nt variants since we already filtered out indels
                ind = mutated_positions.index(position)
                nts_at_allowed_positions[ind] = new_allele

        codons_at_allowed_positions = []
        aas_at_allowed_positions = []
        codon = ""
        for i in range(len(nts_at_allowed_positions)):
            if i % 3 == 2:
                codon += nts_at_allowed_positions[i]
                codons_at_allowed_positions.append(codon)
                aas_at_allowed_positions.append(str(Seq.translate(codon)))
                codon = ""
            else:
                codon += nts_at_allowed_positions[i]
        print(aas_at_allowed_positions)
        var_name = "".join(aas_at_allowed_positions)
        D1135s.append(aas_at_allowed_positions[0])
        S1136s.append(aas_at_allowed_positions[1])
        G1218s.append(aas_at_allowed_positions[2])
        E1219s.append(aas_at_allowed_positions[3])
        R1333s.append(aas_at_allowed_positions[4])
        R1335s.append(aas_at_allowed_positions[5])
        T1337s.append(aas_at_allowed_positions[6])
    else:
        var_name = "NA"
        D1135s.append("NA")
        S1136s.append("NA")
        G1218s.append("NA")
        E1219s.append("NA")
        R1333s.append("NA")
        R1335s.append("NA")
        T1337s.append("NA")
    var_names.append(var_name)

    # add warnings for missense variants where they are not supposed to be
    coding_variant_warning = []
    over_50 = False
    for i in range(len(mismatch_table)):
        row = mismatch_table.iloc[i,]
        if not row["position"] in mutated_positions:
            coding_variant_warning.append(["position: " + str(row["position"]) + " alt_nt: " + str(
                row["alt_nt"]) + " frequency: " + str(row["alt_frequencies"])])
            if row["position"] in recoded_positions and len(row["alt_nt"]) == 1 and row["alt_nt"][0] == \
                    recoded_positions[row["position"]]:
                if any([freq > 0.5 for freq in row["alt_frequencies"]]):
                    over_50 = "intentionally recoded positions only"
            else:
                if any([freq > 0.5 for freq in row["alt_frequencies"]]):
                    over_50 = True
    coding_variant_warnings.append(coding_variant_warning)
    coding_variant_over_50.append(over_50)


df = pd.DataFrame({"sample": samples, "variant_type": var_types, "variant_name": var_names,  "low_reads_warning": low_reads_warnings, "position_of_low_coverage": min_read_positions,
                   "1135": D1135s, "1136": S1136s, "1218": G1218s, "1219": E1219s, "1333": R1333s, "1335": R1335s,
                   "1337": T1337s, "Coding_variant_over_50_percent": coding_variant_over_50, "Coding_variant_warnings": coding_variant_warnings})
df.to_csv("output_variants.csv")












