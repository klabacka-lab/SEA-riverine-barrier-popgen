#!/usr/bin/env python3
"""
Convert interleaved NEXUS alignment files to non-interleaved (sequential) NEXUS.

Usage:
    python make_non_interleaving.py /path/to/nexus_dir

Outputs:
    For each .nex / .nexus file, writes <name>.nex next to it.
"""

from __future__ import annotations

import argparse
import re
from collections import OrderedDict
from pathlib import Path


NEXUS_EXTS = {".nex", ".nexus"}


def strip_inline_comments(s: str) -> str:
    # Simple comment removal for standard [ ... ] comments.
    return re.sub(r"\[.*?\]", "", s)


def parse_matrix_line(line: str):
    """
    Parse one matrix line.

    Returns:
        (label, seq, explicit_label)

    explicit_label=False means this line is treated as an unlabeled continuation
    sequence chunk from an interleaved block.
    """
    line = strip_inline_comments(line).strip()
    if not line:
        return None, None, False

    # Quoted taxon label: 'Taxon name' ACTG...
    m = re.match(r"""^(['"])(.*?)\1\s+(.*)$""", line)
    if m:
        label = m.group(2)
        seq = m.group(3).replace(" ", "").replace("\t", "")
        return label, seq, True

    parts = line.split()
    if len(parts) >= 2:
        # Unquoted label followed by sequence chunk.
        label = parts[0]
        seq = "".join(parts[1:])
        return label, seq, True

    # Single token: usually an unlabeled continuation chunk in later interleaved blocks.
    return None, parts[0], False


def split_blocks(matrix_lines):
    """
    Split matrix lines into blocks separated by blank lines.
    """
    blocks = []
    current = []
    for line in matrix_lines:
        if not line.strip():
            if current:
                blocks.append(current)
                current = []
        else:
            current.append(line)
    if current:
        blocks.append(current)
    return blocks


def parse_interleaved_matrix(matrix_lines):
    """
    Parse interleaved matrix lines and reconstruct full sequences.
    """
    blocks = split_blocks(matrix_lines)
    taxa_order = []
    seqs = OrderedDict()

    for block_index, block in enumerate(blocks):
        parsed = []
        for line in block:
            label, seq, explicit = parse_matrix_line(line)
            if seq is None:
                continue
            parsed.append((label, seq, explicit))

        if not parsed:
            continue

        if block_index == 0:
            # First block should define the taxon order.
            for label, seq, explicit in parsed:
                if not explicit:
                    raise ValueError(
                        "First matrix block looks unlabeled; this script expects a standard interleaved NEXUS file."
                    )
                if label not in taxa_order:
                    taxa_order.append(label)
                seqs.setdefault(label, "")
                seqs[label] += seq
        else:
            # Later blocks: use explicit labels when present, otherwise assume
            # the same taxon order as the first block.
            unlabeled_i = 0
            for label, seq, explicit in parsed:
                if explicit:
                    if label not in seqs:
                        # New label not seen before; keep it in order of appearance.
                        taxa_order.append(label)
                        seqs[label] = ""
                    seqs[label] += seq
                else:
                    if unlabeled_i >= len(taxa_order):
                        raise ValueError(
                            "Found more unlabeled matrix rows than taxa in the first block."
                        )
                    taxon = taxa_order[unlabeled_i]
                    seqs[taxon] += seq
                    unlabeled_i += 1

    return taxa_order, seqs


def quote_label(label: str) -> str:
    """
    Quote a label if needed for NEXUS output.
    """
    if re.fullmatch(r"[A-Za-z0-9._-]+", label):
        return label
    return "'" + label.replace("'", "''") + "'"


def convert_text(text: str) -> str:
    lines = text.splitlines(keepends=True)
    out = []
    i = 0

    while i < len(lines):
        line = lines[i]
        low = line.lower()

        # Normalize FORMAT lines: remove interleave=yes if present.
        if re.match(r"^\s*format\b", line, flags=re.IGNORECASE):
            line = re.sub(r"(?i)\binterleave\s*=\s*yes\b", "interleave=no", line)
            line = re.sub(r"(?i)\binterleave\s*=\s*y\b", "interleave=no", line)
            out.append(line)
            i += 1
            continue

        # Detect MATRIX block.
        if re.match(r"^\s*matrix\b", line, flags=re.IGNORECASE):
            out.append(line)
            i += 1

            matrix_lines = []
            end_semicolon_line = None

            while i < len(lines):
                cur = lines[i]
                cur_nocomment = strip_inline_comments(cur)
                if ";" in cur_nocomment:
                    # End of MATRIX.
                    before, after = cur.split(";", 1)
                    if before.strip():
                        matrix_lines.append(before + "\n")
                    end_semicolon_line = ";" + after
                    i += 1
                    break
                else:
                    matrix_lines.append(cur)
                    i += 1

            taxa_order, seqs = parse_interleaved_matrix(matrix_lines)

            if taxa_order:
                max_label = max(len(quote_label(t)) for t in taxa_order)
                for taxon in taxa_order:
                    label = quote_label(taxon)
                    out.append(f"    {label.ljust(max_label)}  {seqs[taxon]}\n")

            # Close matrix
            if end_semicolon_line is None:
                out.append("    ;\n")
            else:
                out.append("    ;" + end_semicolon_line[1:] if end_semicolon_line.startswith(";") else end_semicolon_line)

            continue

        out.append(line)
        i += 1

    return "".join(out)


def process_file(path: Path) -> Path:
    text = path.read_text(encoding="utf-8")
    converted = convert_text(text)

    out_path = path.with_suffix("")  # drop last suffix
    out_path = out_path.with_name(out_path.name + ".nex")
    out_path.write_text(converted, encoding="utf-8")
    return out_path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", help="Directory containing .nex or .nexus files")
    args = parser.parse_args()

    directory = Path(args.directory).expanduser().resolve()
    if not directory.is_dir():
        raise SystemExit(f"Not a directory: {directory}")

    files = sorted(
        p for p in directory.iterdir()
        if p.is_file() and p.suffix.lower() in NEXUS_EXTS
    )

    if not files:
        print("No .nex or .nexus files found.")
        return

    for f in files:
        try:
            out = process_file(f)
            print(f"Converted: {f.name} -> {out.name}")
        except Exception as e:
            print(f"Skipped {f.name}: {e}")


if __name__ == "__main__":
    main()
