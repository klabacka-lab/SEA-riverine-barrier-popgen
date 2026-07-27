#!/usr/bin/env python3

"""
Create one aligned FASTA file from a CSV containing specimen metadata and GenBank accessions.

Expected CSV columns:
- genus
- specific_epithet
- museum_id
- genbank_id
- river
- side

This script assumes the CSV contains records for a single gene only.

Requirements:
    pip install biopython
    MAFFT installed and available on PATH

NCBI note:
    Set Entrez.email to your email address when fetching from GenBank.

Example use:
    python make_alignment.py \
      --csv specimens.csv \
      --email your.email@example.com \
      --gene-name COI
"""

from __future__ import annotations

import argparse
import csv
import re
import subprocess
import sys
import time
from io import StringIO
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from Bio import Entrez, SeqIO



def sanitize_name(text: str) -> str:
    text = str(text).strip()
    text = re.sub(r"[^\w.-]+", "_", text)
    text = re.sub(r"_+", "_", text)
    return text.strip("_") or "unknown"


def fetch_genbank_fasta(
    accession: str,
    email: str,
    api_key: Optional[str] = None,
    sleep_seconds: float = 0.34,
) -> Tuple[str, str]:
    """
    Fetch a GenBank accession as FASTA text and return (sequence, description).
    """
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    if sleep_seconds > 0:
        time.sleep(sleep_seconds)

    handle = Entrez.efetch(db="nuccore", id=accession, rettype="fasta", retmode="text")
    fasta_text = handle.read()
    handle.close()

    if not fasta_text.strip():
        raise ValueError(f"No FASTA returned for accession {accession}")

    record = SeqIO.read(StringIO(fasta_text), "fasta")
    return str(record.seq), record.description


def run_mafft(input_fasta: Path, output_fasta: Path, mafft_cmd: str = "mafft") -> None:
    """
    Align sequences in input_fasta with MAFFT and write alignment to output_fasta.
    """
    with open(output_fasta, "w") as out_fh:
        result = subprocess.run(
            [mafft_cmd, "--auto", str(input_fasta)],
            stdout=out_fh,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )

    if result.returncode != 0:
        raise RuntimeError(
            f"MAFFT failed for {input_fasta.name}.\n"
            f"stderr:\n{result.stderr}"
        )

def uppercase_fasta_sequences(input_fasta: Path, output_fasta: Path) -> None:
    with open(input_fasta) as inp, open(output_fasta, "w") as out:
        for line in inp:
            if line.startswith(">"):
                out.write(line)
            else:
                out.write(line.upper())


def write_fasta(records: List[Tuple[str, str]], output_path: Path) -> None:
    """
    records: list of (header, sequence)
    """
    with open(output_path, "w") as fh:
        for header, seq in records:
            fh.write(f">{header}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i + 80] + "\n")


def infer_output_basename(rows: List[dict], gene_name: str, genus_col: str, species_col: str, river_col: str) -> str:
    genera = sorted({str(row[genus_col]).strip() for row in rows if str(row[genus_col]).strip()})
    species = sorted({str(row[species_col]).strip() for row in rows if str(row[species_col]).strip()})
    rivers = sorted({str(row[river_col]).strip() for row in rows if str(row[river_col]).strip()})

    genus_part = genera[0] if len(genera) == 1 else "Mixed"
    species_part = species[0] if len(species) == 1 else "spp"
    river_part = rivers[0] if len(rivers) == 1 else "MultipleRivers"

    return sanitize_name(f"{genus_part}_{species_part}_{gene_name}_{river_part}")


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Fetch GenBank accessions from a CSV and create one aligned FASTA file."
    )
    parser.add_argument("--outdir", default="aligned_fastas", help = "output directory for alignments")
    parser.add_argument("--csv", required=True, help="Input CSV file")
    parser.add_argument("--email", required=True, help="NCBI email address required by Entrez")
    parser.add_argument("--api-key", default=None, help="NCBI API key (optional)")
    parser.add_argument("--gene-name", required=True, help="Gene name to use in output file and FASTA headers")
    parser.add_argument("--genbank-column", default="genbank_id", help="GenBank accession column name (default: genbank_id)")
    parser.add_argument("--genus-column", default="genus", help="Genus column name (default: genus)")
    parser.add_argument("--species-column", default="specific_epithet", help="Specific epithet column name (default: specific_epithet)")
    parser.add_argument("--museum-column", default="museum_id", help="Museum ID column name (default: museum_id)")
    parser.add_argument("--river-column", default="river", help="River column name (default: river)")
    parser.add_argument("--side-column", default="side", help="Side column name (default: side)")
    parser.add_argument("--mafft", default="mafft", help="MAFFT command/path (default: mafft)")
    parser.add_argument("--no-align", action="store_true", help="Write unaligned FASTA only; skip MAFFT")
    parser.add_argument("--sleep", type=float, default=0.34, help="Seconds to sleep between GenBank requests (default: 0.34)")
    parser.add_argument(
        "--output-name",
        default=None,
        help="Optional explicit output basename. If omitted, it is inferred as Genus_specificepithet_Gene_River or Genus_spp_Gene_River.",
    )
    args = parser.parse_args()

    csv_path = Path(args.csv)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    with open(csv_path, newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        rows = list(reader)

    if not rows:
        print("CSV is empty.", file=sys.stderr)
        return 1

    required_cols = [
        args.genbank_column,
        args.genus_column,
        args.species_column,
        args.museum_column,
        args.river_column,
        args.side_column,
    ]
    missing = [c for c in required_cols if c not in rows[0]]
    if missing:
        print(f"Missing required columns: {', '.join(missing)}", file=sys.stderr)
        return 1

    output_base = args.output_name or infer_output_basename(
        rows=rows,
        gene_name=args.gene_name,
        genus_col=args.genus_column,
        species_col=args.species_column,
        river_col=args.river_column,
    )

    raw_fasta = outdir / f"{output_base}.unaligned.fasta"
    aln_fasta = outdir / f"{output_base}.fasta"

    accession_cache: Dict[str, Tuple[str, str]] = {}
    fasta_records: List[Tuple[str, str]] = []
    skipped = 0

    for row in rows:
        accession = str(row[args.genbank_column]).strip()
        if not accession:
            skipped += 1
            continue

        try:
            if accession not in accession_cache:
                accession_cache[accession] = fetch_genbank_fasta(
                    accession=accession,
                    email=args.email,
                    api_key=args.api_key,
                    sleep_seconds=args.sleep,
                )

            seq, _desc = accession_cache[accession]

            genus = sanitize_name(row[args.genus_column])
            species = sanitize_name(row[args.species_column])
            museum_id = sanitize_name(row[args.museum_column])
            river = sanitize_name(row[args.river_column])
            side = sanitize_name(row[args.side_column])
            gene_name = sanitize_name(args.gene_name)

            # Header format:
            # >Genus_species_MuseumID_GenbankID_Gene_River_side
            header = "_".join([
                genus,
                species,
                museum_id,
                sanitize_name(accession),
                gene_name,
                river,
                side,
            ])

            fasta_records.append((header, seq))

        except Exception as e:
            skipped += 1
            print(f"[WARN] Skipping {accession}: {e}", file=sys.stderr)

    if not fasta_records:
        print("No sequences were written.", file=sys.stderr)
        return 1

    write_fasta(fasta_records, raw_fasta)
    print(f"Wrote {len(fasta_records)} sequences to {raw_fasta}")

    if args.no_align:
        print("Skipped alignment.")
        return 0

    try:
        temp_aln = outdir / f"{output_base}.aligned.tmp.fasta"
        run_mafft(raw_fasta, temp_aln, mafft_cmd=args.mafft)
        uppercase_fasta_sequences(temp_aln, aln_fasta)
        temp_aln.unlink(missing_ok=True)
        print(f"Aligned FASTA written to {aln_fasta}")
    except FileNotFoundError:
        print(
            f"[ERROR] MAFFT not found: '{args.mafft}'. Install MAFFT or use --no-align.",
            file=sys.stderr,
        )
        return 1
    except Exception as e:
        print(f"[ERROR] Alignment failed: {e}", file=sys.stderr)
        return 1

    if skipped:
        print(f"[INFO] Skipped {skipped} rows", file=sys.stderr)

    print("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
