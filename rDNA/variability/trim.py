from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
import sys

"""
Trim a multiple sequence alignment (MSA) by removing leading and trailing columns
that contain gaps in any sequence.

USAGE
-----
python trim.py <aligned_fasta> <out_fasta>

INPUTS
------
<aligned_fasta> : alignment file (output from MSA tool) in fasta format
<out_fasta>     : output filename for the trimmed output alignment

OUTPUT
------
A trimmed FASTA alignment which starts and ends with fully non-gap columns.
"""

if len(sys.argv) < 3:
    print("Usage: python trim.py <aligned_fasta> <out_fasta>")
    sys.exit(1)

alignment_file = sys.argv[1]
out_file = sys.argv[2]

# load alignment
alignment = AlignIO.read(alignment_file, "fasta")
#print(alignment)

# Find alignment boundaries
alignment_length = alignment.get_alignment_length()
first_non_gap = None
last_non_gap = None

for i in range(alignment_length):
    col = alignment[:, i]
    if all([base != '-' for base in col]):
        first_non_gap = i
        break

for i in range(alignment_length - 1, -1, -1):
    col = alignment[:, i]
    if all([base != '-' for base in col]):
        last_non_gap = i
        break

if first_non_gap is None or last_non_gap is None:
    print("No non-gap columns found!")
    sys.exit(1)

print(f"Trimming alignment ({alignment_length} length): keeping columns {first_non_gap + 1}–{last_non_gap + 1} (1-based) from original alignment.")

# trim sequences
trimmed_records = []
for rec in alignment:
    trimmed_seq = rec.seq[first_non_gap:last_non_gap + 1]
    rec.seq = trimmed_seq
    trimmed_records.append(rec)

# Create new MultipleSeqAlignment object
trimmed_alignment = MultipleSeqAlignment(trimmed_records)

# Save trimmed alignment
AlignIO.write(trimmed_alignment, out_file, "fasta")
print(f"Trimmed alignment saved as {out_file}")
