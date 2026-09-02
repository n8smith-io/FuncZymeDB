#!/usr/bin/env python3
"""
Add NPClassifier annotations from a FuncZymeDB JSON file to a compound TSV.

Usage:
    python code/database/export_compound_classifications.py \
      --json results/database/funczymedb.json \
      --input data/curated/compound_data.tsv \
      --output results/database/compound_data_with_npc.tsv
"""
import argparse
import csv
import json


NPC_FIELDS = (
    ("NPC_class", "class_results"),
    ("NPC_superclass", "superclass_results"),
    ("NPC_pathway", "pathway_results"),
)


def clean(value):
    return str(value or "").strip()


def key(value):
    return clean(value).casefold()


def id_key(value):
    cleaned = clean(value)
    if cleaned.upper().startswith("CHEBI:"):
        cleaned = cleaned.split(":", 1)[1]
    return cleaned.casefold()


def add_lookup(lookup, field, value, compound):
    cleaned = key(value)
    if cleaned:
        lookup[(field, cleaned)] = compound


def build_compound_lookup(compounds):
    lookup = {}
    for compound_name, compound in compounds.items():
        add_lookup(lookup, "name", compound_name, compound)
        add_lookup(lookup, "name", compound.get("canonical_name"), compound)
        for alt_name in compound.get("alternative_names") or []:
            add_lookup(lookup, "name", alt_name, compound)
        add_lookup(lookup, "inchikey", compound.get("inchi_key"), compound)
        add_lookup(lookup, "cid", compound.get("pubchem_id"), compound)
        add_lookup(lookup, "chebi", compound.get("chebi_id"), compound)
    return lookup


def find_compound(row, lookup):
    checks = (
        ("name", key(row.get("COMPOUND_NAME"))),
        ("inchikey", key(row.get("INCHIKEY"))),
        ("cid", id_key(row.get("CID"))),
        ("chebi", id_key(row.get("CHEBI"))),
    )
    for field, value in checks:
        if value and (field, value) in lookup:
            return lookup[(field, value)]
    for name_part in clean(row.get("COMPOUND_NAME")).split(";"):
        value = key(name_part)
        if value and ("name", value) in lookup:
            return lookup[("name", value)]
    return None


def format_npc_values(compound, npc_key):
    if not compound:
        return ""
    npclassifier = compound.get("npclassifier") or {}
    values = []
    seen = set()
    for value in npclassifier.get(npc_key) or []:
        cleaned = clean(value)
        if cleaned and cleaned not in seen:
            values.append(cleaned)
            seen.add(cleaned)
    return " ; ".join(values)


def main():
    parser = argparse.ArgumentParser(
        description="Copy a compound TSV and append NPClassifier category columns."
    )
    parser.add_argument(
        "--json",
        required=True,
        help="FuncZymeDB JSON containing compound npclassifier records.",
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Input compound TSV.",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output compound TSV with NPC columns.",
    )
    args = parser.parse_args()

    with open(args.json) as handle:
        database = json.load(handle)
    lookup = build_compound_lookup(database.get("compounds") or {})

    missing = []
    with open(args.input, newline="") as in_handle, open(
        args.output, "w", newline=""
    ) as out_handle:
        reader = csv.DictReader(in_handle, delimiter="\t")
        fieldnames = list(reader.fieldnames or [])
        for field, _ in NPC_FIELDS:
            if field not in fieldnames:
                fieldnames.append(field)

        writer = csv.DictWriter(
            out_handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()

        for row in reader:
            compound = find_compound(row, lookup)
            if compound is None:
                missing.append(clean(row.get("COMPOUND_NAME")))
            for field, npc_key in NPC_FIELDS:
                row[field] = format_npc_values(compound, npc_key)
            writer.writerow(row)

    print(f"Wrote {args.output}")
    if missing:
        print(f"Rows without JSON compound match: {len(missing)}")
        for compound_name in missing[:20]:
            print(f"  {compound_name}")


if __name__ == "__main__":
    main()
