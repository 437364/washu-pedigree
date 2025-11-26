from Bio import AlignIO
import sys

"""
Compute per-column nucleotide statistics from a multiple FASTA alignment and
identify conserved columns (columns with 100% identity across all sequences).

The script outputs:
1. A TSV file listing non-conserved columns and the counts of A/C/G/T and gaps ('-')
2. Summary statistics printed to stdout

USAGE
-----
python make_statistics.py <aligned_fasta> <stats_file_tsv>
    
INPUTS
------
<aligned_fasta> : alignment file (output from MSA tool) in fasta format
<stats_file_tsv> : path for the output TSV file listing variable columns
   
"""

if len(sys.argv) < 3:
    print("Usage: python make_statistics.py <aligned_fasta> <stats_file_tsv>")
    sys.exit(1)

alignment_file = sys.argv[1]
alignment = AlignIO.read(alignment_file, "fasta")

stats_file = sys.argv[2]

total_length = alignment.get_alignment_length()
conserved_count = 0
conserved_a, conserved_c, conserved_g, conserved_t = 0, 0, 0, 0

with open(stats_file, "w") as out:
    # write header
    out.write("position\ta\tc\tg\tt\tgaps\n")

    for i in range(total_length):
        col = alignment[:, i]
        col = col.lower()

        a, c, g, t = col.count("a"), col.count("c"), col.count("g"), col.count("t")
        gaps = col.count("-")

        col_set = set([base for base in col])
        col_set_arr = list(col_set)
        if len(col_set) == 1:
            conserved_count += 1
            if col_set_arr[0] == "a":
                conserved_a += 1
            elif col_set_arr[0] == "c":
                conserved_c += 1
            elif col_set_arr[0] == "g":
                conserved_g += 1
            elif col_set_arr[0] == "t":
                conserved_t += 1
        else:
            out.write(f"{i + 1}\t{a}\t{c}\t{g}\t{t}\t{gaps}\n")


print(f"Conserved (100% identity) columns: {conserved_count}/{total_length}, {(conserved_count/total_length) * 100}%")
print(f"Conserved cols bases ratio:")
print(f"A: {conserved_a} / {conserved_count}, {(conserved_a / conserved_count) * 100}%")
print(f"C: {conserved_c} / {conserved_count}, {(conserved_c / conserved_count) * 100}%")
print(f"G: {conserved_g} / {conserved_count}, {(conserved_g / conserved_count) * 100}%")
print(f"T: {conserved_t} / {conserved_count}, {(conserved_t / conserved_count) * 100}%")
