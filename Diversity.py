import numpy as np
import pandas as pd
import sys
from collections import Counter
import pysam

# The code takes:
# Reference Sequence Name
# Prefix
# Guide's starting position (count from 0)
# Guide's end position (count from 0)

#
ref_seq_name = sys.argv[1]
prefix = sys.argv[2]
guide_start = sys.argv[3]
guide_end = sys.argv[4]


samfile = pysam.AlignmentFile("%s.bam" % prefix, "rb")
Paired_reads = samfile.count(contig=ref_seq_name) / 2
pysam.faidx("-o", "%s.fa.fai" % ref_seq_name, "%s.fa" % ref_seq_name)
ref = pysam.FastaFile("%s.fa" % ref_seq_name)
seq1 = ref.fetch(reference=ref_seq_name)
print("Paired reads =", Paired_reads)
target_ref = seq1[0:int(guide_start)]
guide_ref = seq1[int(guide_start):int(guide_end)]
threeprime_ref = seq1[int(guide_end):len(seq1)]

all_reads = np.ndarray((samfile.count(), 5), dtype=object)
n = 0
good_reads = np.ndarray((samfile.count(), 1), dtype=object)
i = 0
guide_seq = np.ndarray((samfile.count(), 2), dtype=object)
j = 0
bad_reads = np.ndarray((samfile.count(), 2), dtype=object)
k = 0
indel_reads = np.ndarray((samfile.count(), 1), dtype=object)
m = 0

allreads = samfile.fetch(reference=ref_seq_name, multiple_iterators=True)

for reads in allreads:
    all_reads[n][0] = reads.query_name
    all_reads[n][1] = reads.query_length
    all_reads[n][2] = reads.query_sequence[0:int(guide_start)]
    all_reads[n][3] = reads.query_sequence[int(guide_start):int(guide_end)]
    all_reads[n][4] = reads.query_sequence[int(guide_end):len(seq1)]
    n += 1

reads_matrix = pd.DataFrame(all_reads).dropna()
print("Reads matrix =", len(reads_matrix))
for m in range(len(reads_matrix)):
    if reads_matrix.iloc[m][0] not in indel_reads:
        if reads_matrix.iloc[m][1] != len(seq1):
            indel_reads[m][0] = reads_matrix.iloc[m][0]

print("Target ref =", str(target_ref))
for k in range(len(reads_matrix)):
    if reads_matrix.iloc[k][0] not in indel_reads:
        if reads_matrix.iloc[k][0] not in bad_reads:
            if str(reads_matrix.iloc[k][2]) != str(target_ref):
                print("Target error =", str(reads_matrix.iloc[k][2]))
                error = 0
                for pos in range(len(str(reads_matrix.iloc[k][2]))):
                    base = str(reads_matrix.iloc[k][2])[pos]
                    ref_base = str(target_ref)[pos]
                    if base != ref_base:
                        error += 1
                if error > 6:
                    bad_reads[k][0] = reads_matrix.iloc[k][0]
                    bad_reads[k][1] = error
                elif reads_matrix.iloc[k][4] != threeprime_ref:
                    print("3 prime error =", str(reads_matrix.iloc[k][2]))
                    for pos in range(len(str(reads_matrix.iloc[k][4]))):
                        base = str(reads_matrix.iloc[k][4])[pos]
                        ref_base = str(threeprime_ref)[pos]
                        if base != ref_base:
                            error += 1
                    if error > 6:
                        bad_reads[k][0] = reads_matrix.iloc[k][0]
                        bad_reads[k][1] = error
            elif reads_matrix.iloc[k][4] != threeprime_ref:
                error = 0
                for pos in range(len(str(reads_matrix.iloc[k][4]))):
                    base = str(reads_matrix.iloc[k][4])[pos]
                    ref_base = str(threeprime_ref)[pos]
                    if base != ref_base:
                        error += 1
                if error > 6:
                    bad_reads[k][0] = reads_matrix.iloc[k][0]
                    bad_reads[k][1] = error
                print("Read = ", k, "Errors =", error)

for i in range(len(reads_matrix)):
    if reads_matrix.iloc[i][0] not in indel_reads:
        if reads_matrix.iloc[i][0] not in bad_reads:
            if reads_matrix.iloc[i][0] not in good_reads:
                good_reads[i][0] = reads_matrix.iloc[i][0]
                guide_seq[i][0] = reads_matrix.iloc[i][3]
                guide_seq[i][1] = reads_matrix.iloc[i][0]

indel_matrix = pd.DataFrame(indel_reads).dropna().drop_duplicates()
print("indel matrix =", len(indel_matrix))
error_matrix = pd.DataFrame(bad_reads).dropna()
print("error matrix =", len(error_matrix))
error_type_count = Counter(error_matrix[1])
error_dic = dict(error_type_count)
error_histogram = np.empty((len(error_type_count), 2), dtype=object)
row = 0

for key in error_dic.keys():
    error_histogram[row][0] = key
    error_histogram[row][1] = error_dic[key]
    row += 1

error_Columns = ["# of Mismatches", "Reads"]
error_Event_matrix = pd.DataFrame(error_histogram, columns=error_Columns).sort_values(by="Reads", ascending=False).reset_index(
    drop=True)

guide_matrix = pd.DataFrame(guide_seq).dropna()
seq_type_count = Counter(guide_matrix[0])
total_molecules = len(guide_matrix)
dic = dict(seq_type_count)

freq_histogram = np.empty((len(seq_type_count), 5), dtype=object)
row = 0


logo_histogram = np.zeros((4, (int(int(guide_end)-int(guide_start)))), dtype=int)


for key in dic.keys():
    per = round(float(float(dic[key]) / total_molecules) * 100, 2)
    freq_histogram[row][0] = key
    freq_histogram[row][1] = dic[key]
    freq_histogram[row][2] = per
    MM = 0
    position = str()
    for pos in range(len(str(key))):
        base = str(key)[pos]
        ref_base = str(guide_ref)[pos]
        if base != ref_base:
            MM += 1
            position += str(pos + 1) + ","
        if base == "A":
            logo_histogram[0][pos] += 1
        elif base == "C":
            logo_histogram[1][pos] += 1
        elif base == "G":
            logo_histogram[2][pos] += 1
        elif base == "T":
            logo_histogram[3][pos] += 1
        print(base)
        print(logo_histogram)
    freq_histogram[row][3] = MM
    freq_histogram[row][4] = position[:-1]
    row += 1

Logo_matrix = pd.DataFrame(logo_histogram, index=["A", "C", "G", "T"]).dropna()
Columns = ["Sequence", "Reads", "Percentage", "Mismatch", "Positions"]
Event_matrix = pd.DataFrame(freq_histogram, columns=Columns).sort_values(by="Reads", ascending=False).reset_index(
    drop=True)

mm_type_count = Counter(Event_matrix["Mismatch"])
mmdic = dict(mm_type_count)
unique_molecules = len(Event_matrix)

mm_histogram = np.empty((len(mm_type_count), 3), dtype=object)
mmrow = 0

for mmkey in mmdic.keys():
    mmper = round(float(float(mmdic[mmkey]) / unique_molecules) * 100, 2)
    mm_histogram[mmrow][0] = mmkey
    mm_histogram[mmrow][1] = mmdic[mmkey]
    mm_histogram[mmrow][2] = mmper
    mmrow += 1

mm_Columns = ["Number of Mismatches", "Molecules", "Percentage"]
mm_matrix = pd.DataFrame(mm_histogram, columns=mm_Columns).sort_values(by="Molecules", ascending=False).reset_index(
    drop=True)

Stats_data = np.empty((5, 3), dtype=object)
Stats_data[0][0] = "Initial paired reads"
Stats_data[0][1] = Paired_reads
Stats_data[0][2] = round(float(float(Paired_reads) / Paired_reads) * 100, 2)
Stats_data[1][0] = "Reads with Indels"
Stats_data[1][1] = len(indel_matrix)
Stats_data[1][2] = round(float(float(len(indel_matrix)) / Paired_reads) * 100, 2)
Stats_data[2][0] = "Reads with mismatches"
Stats_data[2][1] = len(error_matrix)
Stats_data[2][2] = round(float(float(len(error_matrix)) / Paired_reads) * 100, 2)
Stats_data[3][0] = "Reads analyzed"
Stats_data[3][1] = len(guide_matrix)
Stats_data[3][2] = round(float(float(len(guide_matrix)) / Paired_reads) * 100, 2)
Stats_data[4][0] = "Different guide sequences"
Stats_data[4][1] = len(seq_type_count)
Stats_data[4][2] = round(float(float(len(seq_type_count)) / len(guide_matrix)) * 100, 2)

Stat_Columns = ["Parameter", "Total", "Percentage"]
Stats_data = pd.DataFrame(Stats_data, columns=Stat_Columns).reset_index(
    drop=True)

with pd.ExcelWriter("%s_Diversity.xlsx" % prefix) as writer:
    Event_matrix.to_excel(writer, sheet_name="Diversity")
    Stats_data.to_excel(writer, sheet_name="Stats")
    mm_matrix.to_excel(writer, sheet_name="Guide_MM_histogram")
    error_Event_matrix.to_excel(writer, sheet_name="error_MM_histogram")
    Logo_matrix.to_excel(writer, sheet_name="Logo_matrix")