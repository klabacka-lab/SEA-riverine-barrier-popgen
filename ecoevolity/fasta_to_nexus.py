from Bio import AlignIO
import os

# Text file containing one FASTA filename per line
list_file = "alignment_list.txt"

with open(list_file) as f:
    fasta_files = [line.strip() for line in f if line.strip()]

for fasta_file in fasta_files:
    try:
        # Read the alignment
        alignment = AlignIO.read(fasta_file, "fasta")

        # Specify molecule type
        for record in alignment:
            record.annotations["molecule_type"] = "DNA"

        # Output filename
        nexus_file = os.path.splitext(fasta_file)[0] + ".nex"

        # Write NEXUS file
        AlignIO.write(alignment, nexus_file, "nexus")

        print(f"Converted: {fasta_file} -> {nexus_file}")

    except Exception as e:
        print(f"Error converting {fasta_file}: {e}")
