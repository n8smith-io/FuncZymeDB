#!/usr/bin/env python3
"""
merge_species_lineages.py

Append any species from a one‐per‐line text file into the species_lineages.tsv
if they don’t already appear there, leaving the other columns blank for manual curation.
"""

import argparse
import csv
import os

def parse_args():
    p = argparse.ArgumentParser(
        description="Merge unique species into species_lineages.tsv"
    )
    p.add_argument(
        "-t", "--tsv", required=True,
        help="Path to existing species_lineages.tsv"
    )
    p.add_argument(
        "-u", "--unique", required=True,
        help="Path to text file of unique species (one per line)"
    )
    p.add_argument(
        "-o", "--output", default=None,
        help="Output TSV path (default: same basename as --tsv + '.updated.tsv')"
    )
    return p.parse_args()

def main():
    args = parse_args()

    # determine output path
    if args.output:
        out_tsv = args.output
    else:
        base, _ = os.path.splitext(args.tsv)
        out_tsv = f"{base}.updated.tsv"

    # load existing TSV
    with open(args.tsv, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        fieldnames = reader.fieldnames
        rows = list(reader)

    existing = { row["SPECIES"] for row in rows }

    # read unique list and find missing
    missing = []
    with open(args.unique) as f:
        for line in f:
            name = line.strip()
            if not name or name in existing:
                continue
            missing.append(name)

    # append blank‐filled rows for each missing species
    for name in missing:
        blank_row = { col: "" for col in fieldnames }
        blank_row["SPECIES"] = name
        rows.append(blank_row)

    # write out the merged TSV
    with open(out_tsv, "w", newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)

    print(f"Appended {len(missing)} new species. Wrote output to {out_tsv}")

if __name__ == "__main__":
    main()
