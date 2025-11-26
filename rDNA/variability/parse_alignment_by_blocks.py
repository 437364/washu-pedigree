from Bio import AlignIO
from collections import Counter
import sys

"""
Parse a multiple sequence alignment (MSA) in FASTA format and segment the
alignment into *conserved* and *variable* blocks based on a conservation
threshold.

A column is considered conserved if the most common nucleotide appears with
frequency >= THRESHOLD (gaps '-' are ignored). Variable blocks also report the
counts of nucleotides within the block.

NOTES
--------
- Coordinates are 0-based and inclusive.
- Conserved blocks have no base counts (only 3 columns are written).
- Only variable blocks include counts for A/C/G/T.
- Gaps ("-") are ignored when determining the dominant base.

USAGE
--------
python parse_alignment_by_blocks.py <alignment.fasta> <output.tsv> <threshold>
    
ARGUMENTS
--------
<aligned_fasta>  : alignment file (output from MSA tool) in fasta format
<output.tsv> : Path to the output TSV file.
<conservation_threshold> : Minimum fraction required for a column to be considered conserved.

"""

if len(sys.argv) < 3:
    print("Usage: python parse_alignment_by_blocks.py <alignment.fasta> <output.tsv> <conservation_threshold>")
    sys.exit(1)

alignment_file = sys.argv[1]
alignment = AlignIO.read(alignment_file, "fasta")
seqs_cnt = len(alignment)
aln_len = alignment.get_alignment_length()

output_filename = sys.argv[2]
output_file = open(output_filename, "w")
output_file.write(f"start\tend\ttype\tA_count\tC_count\tG_count\tT_count\n")

try:
    THRESHOLD = int(sys.argv[3])
except TypeError:
    print("Threshold must be an integer.")
    sys.exit(1)

def col_is_conserved(col_data):
    global THRESHOLD
    counts = Counter(col_data)

    # ignore gaps
    if "-" in counts:
        del counts["-"]

    if not counts:
        return False, ""  # no valid bases

    primary_base, primary_count = counts.most_common(1)[0]

    freq = primary_count / len(col_data)
    return freq >= THRESHOLD, primary_base


def get_block_counts(data):
    counts = {
        "a": 0,
        "c": 0,
        "g": 0,
        "t": 0,
        "-": 0
    }

    for base in data:
        counts[base] += 1

    return counts["a"], counts["c"], counts["g"], counts["t"]


curr_status = None  # 0 = conserved, 1 = variable
primary_bases = []
start, end = 0, 0
for i in range(aln_len):
    col = alignment[:, i]
    is_conserved, primary_base = col_is_conserved(col)

    if curr_status is None:
        if is_conserved:
            curr_status = 0
        else:
            curr_status = 1

    elif curr_status == 0:
        if is_conserved: # was conserved, still is conserved
            end += 1
        else: # was conserved, now isn't
            output_file.write(f"{start}\t{end}\tconserved\n")

            end += 1
            start = end

            primary_bases = [primary_base]
            curr_status = 1

    elif curr_status == 1:
        if not is_conserved:  # wasn't conserved and still isn't
            primary_bases.append(primary_base)
            end += 1
        else: # wasn't conserved, now is
            a, c, g, t = get_block_counts(primary_bases)
            output_file.write(f"{start}\t{end}\tvariable\t{a}\t{c}\t{g}\t{t}\n")

            end += 1
            start = end

            curr_status = 0

# final block
if curr_status == 0:
    output_file.write(f"{start}\t{end}\tconserved\n")
else:
    a, c, g, t = get_block_counts(primary_bases)
    output_file.write(f"{start}\t{end}\tvariable\t{a}\t{c}\t{g}\t{t}\n")


output_file.close()