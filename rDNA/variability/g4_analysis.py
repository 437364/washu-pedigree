from Bio import AlignIO
import sys
import pandas as pd

"""
This script maps G4 motif coordinates (from a GFF file) onto a multiple
sequence alignment (MSA) and count conservation score for each G4.

IMPORTANT
-----
aligned_fasta should be the full, non-trimmed alignment so that positional offsets remain correct.

NOTES
-----
- Coordinates in the GFF file are assumed to be 1-based and inclusive.

USAGE
-----
    python g4_analysis.py <aligned_fasta> <g4_data_gff> <output_tsv>
    
ARGUMENTS
---------
    <aligned_fasta>  : Alignment file (output from MSA tool) in fasta format.
    <g4_data_gff> : GFF3 file containing predicted/found G4 features. Required fields: seqid, start, end, score, attributes.
    <output_tsv> : Output TSV file summarizing G4 features with added conservation scores to each G4.
"""

def conserved_col(idx):
    col = alignment[:, idx]
    return len(set([base for base in col])) == 1


def count_idx(alignment, idx, seqid):
    """Walks along the aligned sequence to translate unaligned coordinates into aligned positions."""
    aligned_seq = next((rec for rec in alignment if rec.id == seqid), None)
    res = 0
    while idx > 0:
        if aligned_seq[res] != "-":
            idx -= 1
        res += 1

    return res, aligned_seq[res]

if len(sys.argv) < 4:
    print("Usage: python g4_analysis.py <aligned_fasta> (not trimmed!) <g4_data_gff> <output_tsv>")
    sys.exit(1)

alignment_file = sys.argv[1]
alignment = AlignIO.read(alignment_file, "fasta")

g4_data = sys.argv[2]
output_file = sys.argv[3]

total_length = alignment.get_alignment_length()
conserved_count = 0

cols = [
    "seqid", "source", "type", "start", "end",
    "score", "strand", "phase", "attributes"
]

gff = pd.read_csv(g4_data, sep="\t", comment="#", names=cols)
with open(output_file, "w") as out, open(g4_data, "r") as g4_data:
    # write header
    out.write("seqid\tstart\tend\tscore\tattributes\tsequence_from_alignment\tconserved(%)\n")

    # process each alignment-mapped position
    for _, row in gff.iterrows():
        start, end = row["start"], row["end"]
        length = end - start + 1
        conserved, cnt = 0, 0
        seq_from_alignment = ""
        for idx in range(length):
            curr_idx, curr_nucleotide = count_idx(alignment, start + idx - 1, row["seqid"])
            cnt += 1
            conserved += conserved_col(curr_idx)
            seq_from_alignment += curr_nucleotide.upper()

        out.write(f"{row['seqid']}\t{start}\t{end}\t{row['score']}\t{row['attributes']}\t{seq_from_alignment}\t{conserved / cnt}\n")
        print(start, end, "done")

print(f"Done. Output in {output_file}.")
