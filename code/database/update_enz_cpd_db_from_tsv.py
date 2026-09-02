#!/usr/bin/env python3
"""
Incrementally update an existing enzyme/compound database JSON from updated TSVs.

This script is intended for small, approval-driven repairs where rebuilding the
entire database is unnecessary. It:

1. Reads an existing DB JSON with `enzymes`, `compounds`, and `reactions`.
2. Reads updated compound/activity TSVs.
3. Resolves name-only activity entries against the updated compound TSV.
4. Previews compound-centered change packages and asks for per-package approval.
5. Applies approved updates, optionally fetching NPClassifier annotations only
   for compounds that newly gained usable SMILES and still lack that data.
6. Appends a detailed log entry into the DB JSON.
"""

from __future__ import annotations

import argparse
import copy
import csv
import json
import sys
import tempfile
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Set, Tuple

from build_enz_cpd_seq_database import (
    NA,
    NPClassifierClient,
    get_sanitized_value,
    normalize_smiles,
    sanitize_npclassifier,
    sanitize_string,
    select_smiles,
)


ROLE_COLUMN_MAP = {
    "SUBSTRATE": "substrates",
    "DONOR": "donors",
    "ACCEPTOR": "acceptors",
    "PRODUCT": "products",
}

COMPOUND_DIFF_FIELDS = [
    "canonical_name",
    "alternative_names",
    "inchi_key",
    "pubchem_id",
    "chebi_id",
    "smiles",
    "aromatic",
    "aliphatic",
    "is_donor",
    "is_acceptor",
    "is_product",
    "compound_notes",
    "verified_by",
    "row_index",
    "npclassifier",
]

REACTION_DIFF_FIELDS = [
    "substrates",
    "donors",
    "acceptors",
    "products",
    "reaction_specific_id",
]

ENZYME_DIFF_FIELDS = [
    "associated_compounds",
]


def load_json(path: Path) -> Dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def write_json_atomic(path: Path, payload: Dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        "w",
        encoding="utf-8",
        dir=str(path.parent),
        delete=False,
    ) as handle:
        json.dump(payload, handle, indent=2)
        handle.write("\n")
        temp_name = handle.name
    Path(temp_name).replace(path)


def ordered_unique(values: Iterable[str]) -> List[str]:
    seen: Set[str] = set()
    out: List[str] = []
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        out.append(value)
    return out


def sorted_unique(values: Iterable[str]) -> List[str]:
    return sorted({value for value in values if value is not None})


def normalize_name(name: str) -> str:
    return sanitize_string(name, default="").lower()


def normalize_chebi(value: str) -> str:
    cleaned = sanitize_string(value)
    if cleaned != NA and cleaned.upper().startswith("CHEBI:"):
        cleaned = cleaned.split("CHEBI:", 1)[1]
    return sanitize_string(cleaned)


def normalize_pubchem(value: str) -> str:
    cleaned = sanitize_string(value)
    if cleaned != NA and cleaned.upper().startswith("CID:"):
        cleaned = cleaned.split("CID:", 1)[1]
    return sanitize_string(cleaned)


def default_compound_record(canonical_name: str) -> Dict[str, Any]:
    return {
        "canonical_name": canonical_name,
        "alternative_names": [canonical_name],
        "inchi_key": NA,
        "pubchem_id": NA,
        "chebi_id": NA,
        "smiles": NA,
        "aromatic": False,
        "aliphatic": False,
        "is_donor": False,
        "is_acceptor": False,
        "is_product": False,
        "compound_notes": NA,
        "verified_by": NA,
        "row_index": [],
        "npclassifier": {},
    }


def parse_compound_templates(path: Path) -> Tuple[Dict[str, Dict[str, Any]], Dict[str, str]]:
    templates: Dict[str, Dict[str, Any]] = {}
    alias_lookup: Dict[str, str] = {}
    with path.open("r", newline="", encoding="utf-8-sig", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            raw_name_field = sanitize_string(row.get("COMPOUND_NAME", ""), default="")
            names = [
                sanitize_string(part, default="")
                for part in raw_name_field.split(";")
            ]
            names = [name for name in names if name]
            if not names:
                continue
            canonical = names[0]
            smiles_value, _smiles_verified = select_smiles(
                row.get("SMILES", NA),
                row.get("VERIFIED_SMILES", NA),
            )
            compound_type = get_sanitized_value(row, "COMPOUND_TYPE")
            template = {
                "canonical_name": canonical,
                "alternative_names": sorted_unique(names),
                "inchi_key": get_sanitized_value(row, "INCHIKEY"),
                "pubchem_id": normalize_pubchem(row.get("CID", NA)),
                "chebi_id": normalize_chebi(row.get("CHEBI", NA)),
                "smiles": smiles_value,
                "aromatic": get_sanitized_value(row, "AROMATIC", default="0") == "1",
                "aliphatic": get_sanitized_value(row, "ALIPHATIC", default="0") == "1",
                "is_donor": compound_type.strip().lower() == "donor",
                "is_acceptor": compound_type.strip().lower() == "acceptor",
                "is_product": compound_type.strip().lower() == "product",
                "compound_notes": get_sanitized_value(row, "COMPOUND_NOTES"),
                "verified_by": get_sanitized_value(row, "VERIFIED_BY"),
                "npclassifier": {},
            }
            templates[canonical] = template
            for name in template["alternative_names"]:
                alias_lookup[normalize_name(name)] = canonical
    return templates, alias_lookup


def build_db_compound_alias_lookup(compounds: Dict[str, Dict[str, Any]]) -> Dict[str, str]:
    alias_lookup: Dict[str, str] = {}
    for key, record in compounds.items():
        names = set(record.get("alternative_names") or [])
        names.add(record.get("canonical_name") or key)
        names.add(key)
        for name in names:
            normalized = normalize_name(str(name))
            if normalized:
                alias_lookup[normalized] = key
    return alias_lookup


def parse_activity_tokens(raw_value: str) -> List[str]:
    cleaned = sanitize_string(raw_value, default="")
    if not cleaned:
        return []
    if cleaned.upper() == NA:
        return [NA]
    tokens = [
        sanitize_string(part, default="")
        for part in cleaned.split(";")
    ]
    return [token for token in tokens if token]


def template_from_annotated_token(token: str) -> Optional[Dict[str, Any]]:
    parts = [sanitize_string(part, default="") for part in token.split("|")]
    if len(parts) != 4:
        return None
    name = sanitize_string(parts[0])
    if name == NA:
        return None
    return {
        "canonical_name": name,
        "alternative_names": [name],
        "inchi_key": sanitize_string(parts[2]),
        "pubchem_id": normalize_pubchem(parts[1]),
        "chebi_id": normalize_chebi(parts[3]),
        "smiles": NA,
        "aromatic": False,
        "aliphatic": False,
        "is_donor": False,
        "is_acceptor": False,
        "is_product": False,
        "compound_notes": NA,
        "verified_by": NA,
        "npclassifier": {},
    }


def merge_template_with_annotation(
    template: Dict[str, Any],
    annotated: Optional[Dict[str, Any]],
) -> Dict[str, Any]:
    merged = copy.deepcopy(template)
    if not annotated:
        return merged
    merged["alternative_names"] = sorted_unique(
        list(merged.get("alternative_names") or []) + list(annotated.get("alternative_names") or [])
    )
    for field in ("inchi_key", "pubchem_id", "chebi_id", "smiles"):
        if sanitize_string(str(merged.get(field, NA))) == NA:
            merged[field] = annotated.get(field, NA)
    return merged


def resolve_token_template(
    token: str,
    compound_templates: Dict[str, Dict[str, Any]],
    update_alias_lookup: Dict[str, str],
    db_alias_lookup: Dict[str, str],
    db_compounds: Dict[str, Dict[str, Any]],
) -> Tuple[Optional[str], Optional[Dict[str, Any]], Optional[str]]:
    if token == NA:
        return NA, None, None

    annotated_template = template_from_annotated_token(token)
    candidate_name = None
    if annotated_template:
        candidate_name = normalize_name(annotated_template["canonical_name"])
    else:
        candidate_name = normalize_name(token)

    if candidate_name in update_alias_lookup:
        canonical = update_alias_lookup[candidate_name]
        template = merge_template_with_annotation(
            compound_templates[canonical],
            annotated_template,
        )
        return canonical, template, None

    if candidate_name in db_alias_lookup:
        db_key = db_alias_lookup[candidate_name]
        canonical = db_compounds[db_key].get("canonical_name") or db_key
        return canonical, None, None

    if annotated_template:
        canonical = annotated_template["canonical_name"]
        return canonical, annotated_template, None

    return None, None, token


def diff_fields(before: Optional[Dict[str, Any]], after: Optional[Dict[str, Any]], fields: Sequence[str]) -> Dict[str, Dict[str, Any]]:
    before_obj = before or {}
    after_obj = after or {}
    out: Dict[str, Dict[str, Any]] = {}
    for field in fields:
        before_value = before_obj.get(field)
        after_value = after_obj.get(field)
        if before_value != after_value:
            out[field] = {"before": before_value, "after": after_value}
    return out


def append_missing_preserve_order(base: Sequence[str], additions: Iterable[str]) -> List[str]:
    out = list(base)
    seen = set(out)
    for value in additions:
        if value in seen:
            continue
        out.append(value)
        seen.add(value)
    return out


def compound_before_after(
    db: Dict[str, Any],
    canonical: str,
    template: Dict[str, Any],
    extra_row_indices: Optional[Iterable[str]] = None,
) -> Tuple[Optional[str], Optional[Dict[str, Any]], Dict[str, Any]]:
    compounds = db["compounds"]
    db_alias_lookup = build_db_compound_alias_lookup(compounds)
    matched_keys = ordered_unique(
        db_alias_lookup[name_key]
        for name in [canonical] + list(template.get("alternative_names") or [])
        for name_key in [normalize_name(name)]
        if name_key and name_key in db_alias_lookup
    )

    existing_key: Optional[str] = None
    if canonical in compounds:
        existing_key = canonical
    elif matched_keys:
        existing_key = matched_keys[0]

    before = copy.deepcopy(compounds[existing_key]) if existing_key else None
    base = copy.deepcopy(before) if before else default_compound_record(canonical)
    after = copy.deepcopy(base)
    after["canonical_name"] = canonical
    after["alternative_names"] = append_missing_preserve_order(
        list(base.get("alternative_names") or []),
        list(template.get("alternative_names") or []) + [canonical],
    )
    for field in ("inchi_key", "pubchem_id", "chebi_id", "smiles", "compound_notes", "verified_by"):
        candidate = template.get(field, NA)
        if sanitize_string(str(candidate), default=NA) != NA:
            after[field] = candidate
    for field in ("aromatic", "aliphatic", "is_donor", "is_acceptor", "is_product"):
        after[field] = bool(base.get(field)) or bool(template.get(field))
    after["row_index"] = append_missing_preserve_order(
        list(base.get("row_index") or []),
        list(extra_row_indices or []),
    )
    after["npclassifier"] = sanitize_npclassifier(base.get("npclassifier"))
    return existing_key, before, after


def token_matches_candidate(token: str, candidate_names: Set[str]) -> bool:
    if token == NA:
        return False
    annotated = template_from_annotated_token(token)
    if annotated:
        token_name = normalize_name(annotated["canonical_name"])
    else:
        token_name = normalize_name(token)
    return token_name in candidate_names


def canonicalize_existing_values(
    values: Sequence[str],
    canonical: str,
    replacement_names: Set[str],
) -> List[str]:
    out: List[str] = []
    for value in values:
        normalized = normalize_name(str(value))
        if normalized in replacement_names:
            out.append(canonical)
        else:
            out.append(value)
    return ordered_unique(out)


def build_reaction_change(
    db: Dict[str, Any],
    activity_row: Dict[str, str],
    canonical: str,
    candidate_names: Set[str],
    replacement_names: Set[str],
) -> Tuple[Optional[Dict[str, Any]], Set[str], List[str]]:
    row_index = get_sanitized_value(activity_row, "INDEX")
    current = db["reactions"].get(row_index)
    if current is None:
        return None, set(), []

    unresolved_tokens: List[str] = []
    referenced_canonicals: Set[str] = set()

    proposed = copy.deepcopy(current)
    changed = False

    for tsv_col, json_col in ROLE_COLUMN_MAP.items():
        tokens = parse_activity_tokens(activity_row.get(tsv_col, NA))
        current_values = list(current.get(json_col) or [])
        if tokens == [NA] or not tokens:
            proposed[json_col] = canonicalize_existing_values(current_values, canonical, replacement_names)
            if proposed[json_col] != current_values:
                changed = True
            continue

        proposed_values = canonicalize_existing_values(current_values, canonical, replacement_names)
        matched_token = False
        for token in tokens:
            if token_matches_candidate(token, candidate_names):
                matched_token = True
                referenced_canonicals.add(canonical)
        if matched_token and canonical not in proposed_values:
            proposed_values.append(canonical)

        proposed_values = ordered_unique(proposed_values)
        if proposed_values != current_values:
            proposed[json_col] = proposed_values
            changed = True
        else:
            proposed[json_col] = current_values

    if changed:
        participant_set: Set[str] = set()
        for json_col in ROLE_COLUMN_MAP.values():
            for value in proposed.get(json_col) or []:
                if value != NA:
                    participant_set.add(value)
        sorted_names = sorted(participant_set)
        proposed["reaction_specific_id"] = "_+_".join(sorted_names) if sorted_names else NA
    else:
        proposed["reaction_specific_id"] = current.get("reaction_specific_id", NA)

    diff = diff_fields(current, proposed, REACTION_DIFF_FIELDS)
    if not diff:
        return None, referenced_canonicals, unresolved_tokens

    return {
        "row_index": row_index,
        "enzyme_id": current.get("enzyme_id"),
        "before": current,
        "after": proposed,
        "diff": diff,
    }, referenced_canonicals, unresolved_tokens


def recompute_associated_compounds(db: Dict[str, Any], enzyme_id: str) -> List[str]:
    compounds: List[str] = []
    for reaction_id in db["enzymes"][enzyme_id].get("reactions") or []:
        reaction = db["reactions"].get(reaction_id)
        if not reaction:
            continue
        for json_col in ROLE_COLUMN_MAP.values():
            compounds.extend(reaction.get(json_col) or [])
    return ordered_unique(name for name in compounds if name)


def build_enzyme_change(
    db: Dict[str, Any],
    enzyme_id: str,
    canonical: str,
    replacement_names: Set[str],
) -> Optional[Dict[str, Any]]:
    current = db["enzymes"].get(enzyme_id)
    if current is None:
        return None
    proposed = copy.deepcopy(current)
    current_values = list(current.get("associated_compounds") or [])
    proposed_values = canonicalize_existing_values(current_values, canonical, replacement_names)
    for reaction_id in current.get("reactions") or []:
        reaction = db["reactions"].get(reaction_id)
        if not reaction:
            continue
        for json_col in ROLE_COLUMN_MAP.values():
            proposed_values = append_missing_preserve_order(proposed_values, reaction.get(json_col) or [])
    proposed["associated_compounds"] = ordered_unique(proposed_values)
    diff = diff_fields(current, proposed, ENZYME_DIFF_FIELDS)
    if not diff:
        return None
    return {
        "enzyme_id": enzyme_id,
        "before": current,
        "after": proposed,
        "diff": diff,
    }


def format_value(value: Any) -> str:
    if isinstance(value, str):
        return value
    return json.dumps(value, sort_keys=True)


def print_compound_preview(
    canonical: str,
    compound_change: Dict[str, Any],
    reaction_changes: List[Dict[str, Any]],
    enzyme_changes: List[Dict[str, Any]],
    npclassifier_needed: bool,
    unresolved_tokens: List[str],
) -> None:
    action = "create" if compound_change["before"] is None else "update"
    rename_note = ""
    if compound_change["existing_key"] and compound_change["existing_key"] != canonical:
        rename_note = f" (matched existing key: {compound_change['existing_key']})"
    print()
    print("=" * 88)
    print(f"Compound package: {canonical} [{action}]{rename_note}")
    print("-" * 88)
    if compound_change["diff"]:
        print("Compound changes:")
        for field, payload in compound_change["diff"].items():
            print(f"  - {field}: {format_value(payload['before'])} -> {format_value(payload['after'])}")
    else:
        print("Compound changes: none")

    if reaction_changes:
        print("Reaction changes:")
        for change in reaction_changes:
            print(f"  - {change['row_index']} ({change['enzyme_id']})")
            for field, payload in change["diff"].items():
                print(f"      {field}: {format_value(payload['before'])} -> {format_value(payload['after'])}")
    else:
        print("Reaction changes: none")

    if enzyme_changes:
        print("Enzyme changes:")
        for change in enzyme_changes:
            print(f"  - {change['enzyme_id']}")
            for field, payload in change["diff"].items():
                print(f"      {field}: {format_value(payload['before'])} -> {format_value(payload['after'])}")
    else:
        print("Enzyme changes: none")

    if npclassifier_needed:
        print("NPClassifier: will fetch after approval")
    else:
        print("NPClassifier: no fetch needed")

    if unresolved_tokens:
        print(f"Unresolved activity tokens retained as-is: {json.dumps(sorted_unique(unresolved_tokens))}")


def prompt_for_approval(canonical: str) -> str:
    while True:
        choice = input(f"Apply package for {canonical}? [y]es/[n]o/[a]ll/[q]uit: ").strip().lower()
        if choice in {"y", "n", "a", "q"}:
            return choice
        print("Please enter y, n, a, or q.")


def apply_compound_change(
    db: Dict[str, Any],
    canonical: str,
    compound_change: Dict[str, Any],
) -> None:
    compounds = db["compounds"]
    existing_key = compound_change["existing_key"]
    after = copy.deepcopy(compound_change["after"])
    if existing_key and existing_key != canonical and existing_key in compounds:
        del compounds[existing_key]
    compounds[canonical] = after


def apply_reaction_change(
    db: Dict[str, Any],
    reaction_change: Dict[str, Any],
) -> None:
    db["reactions"][reaction_change["row_index"]] = copy.deepcopy(reaction_change["after"])


def apply_enzyme_change(
    db: Dict[str, Any],
    enzyme_change: Dict[str, Any],
) -> None:
    db["enzymes"][enzyme_change["enzyme_id"]] = copy.deepcopy(enzyme_change["after"])


def needs_npclassifier_fetch(compound_after: Dict[str, Any]) -> bool:
    smiles = normalize_smiles(compound_after.get("smiles", NA))
    npclassifier = sanitize_npclassifier(compound_after.get("npclassifier"))
    return smiles != NA and not npclassifier


def find_candidate_compounds(
    db: Dict[str, Any],
    activity_path: Path,
    compound_templates: Dict[str, Dict[str, Any]],
    update_alias_lookup: Dict[str, str],
) -> Dict[str, List[Dict[str, str]]]:
    candidates: Dict[str, List[Dict[str, str]]] = defaultdict(list)
    with activity_path.open("r", newline="", encoding="utf-8-sig", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            row_index = get_sanitized_value(row, "INDEX")
            if row_index not in db.get("reactions", {}):
                continue
            row_canonicals: Set[str] = set()
            for tsv_col in ROLE_COLUMN_MAP.keys():
                for token in parse_activity_tokens(row.get(tsv_col, NA)):
                    canonical, _template, _unresolved = resolve_token_template(
                        token,
                        compound_templates,
                        update_alias_lookup,
                        build_db_compound_alias_lookup(db["compounds"]),
                        db["compounds"],
                    )
                    if canonical and canonical in compound_templates:
                        row_canonicals.add(canonical)
            for canonical in row_canonicals:
                candidates[canonical].append(row)

    db_alias_lookup = build_db_compound_alias_lookup(db["compounds"])
    for canonical, template in compound_templates.items():
        if canonical in candidates:
            continue
        for name in [canonical] + list(template.get("alternative_names") or []):
            if normalize_name(name) in db_alias_lookup:
                candidates[canonical] = []
                break
    return dict(candidates)


def collect_changed_enzymes(reaction_changes: Sequence[Dict[str, Any]]) -> List[str]:
    return ordered_unique(
        change["enzyme_id"]
        for change in reaction_changes
        if change.get("enzyme_id")
    )


def append_log_entry(
    db: Dict[str, Any],
    args: argparse.Namespace,
    approved_logs: List[Dict[str, Any]],
    skipped: List[str],
) -> None:
    log_entry = {
        "time": datetime.now().isoformat(),
        "person": args.person,
        "input": [
            ["db_json", str(args.db_json)],
            ["compound_tsv", str(args.compound_tsv)],
            ["activity_tsv", str(args.activity_tsv)],
        ],
        "output": [
            ["output_db_json", str(args.output_json)],
        ],
        "approved_packages": approved_logs,
        "skipped_packages": skipped,
    }
    log_section = db.setdefault("log", {})
    prefix = "tsv_incremental_update_"
    existing_indices: List[int] = []
    for key in log_section.keys():
        if not key.startswith(prefix):
            continue
        suffix = key[len(prefix):]
        if suffix.isdigit():
            existing_indices.append(int(suffix))
    next_index = (max(existing_indices) + 1) if existing_indices else 0
    log_section[f"{prefix}{next_index}"] = log_entry


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Incrementally update an existing enzyme/compound DB JSON from updated TSVs."
    )
    parser.add_argument("--db_json", type=Path, required=True, help="Existing DB JSON to update.")
    parser.add_argument(
        "--output_json",
        type=Path,
        default=None,
        help="Output DB JSON path. Defaults to in-place update of --db_json.",
    )
    parser.add_argument("--compound_tsv", type=Path, required=True, help="Updated compound TSV.")
    parser.add_argument("--activity_tsv", type=Path, required=True, help="Updated activity TSV.")
    parser.add_argument("--person", required=True, help="Initials or name for DB logging.")
    parser.add_argument(
        "--compound_name",
        action="append",
        default=[],
        help=(
            "Restrict updates to one compound canonical name or alias from the update TSV. "
            "Repeat this flag to target multiple compounds."
        ),
    )
    parser.add_argument("--yes", action="store_true", help="Apply all detected packages without prompting.")
    parser.add_argument(
        "--no_npclassifier",
        action="store_true",
        help="Do not call NPClassifier, even when new SMILES are added.",
    )
    return parser


def main() -> None:
    parser = build_arg_parser()
    args = parser.parse_args()
    args.output_json = args.output_json or args.db_json

    db = load_json(args.db_json)
    if not {"enzymes", "compounds", "reactions"}.issubset(db.keys()):
        sys.exit("Input DB JSON must contain top-level keys: enzymes, compounds, reactions.")

    compound_templates, update_alias_lookup = parse_compound_templates(args.compound_tsv)
    candidate_rows = find_candidate_compounds(
        db,
        args.activity_tsv,
        compound_templates,
        update_alias_lookup,
    )

    if args.compound_name:
        requested_canonicals: List[str] = []
        missing_names: List[str] = []
        for raw_name in args.compound_name:
            normalized = normalize_name(raw_name)
            canonical = update_alias_lookup.get(normalized)
            if canonical is None and raw_name in compound_templates:
                canonical = raw_name
            if canonical is None:
                missing_names.append(raw_name)
                continue
            requested_canonicals.append(canonical)
        if missing_names:
            sys.exit(
                "The following --compound_name values were not found in the update TSV: "
                + ", ".join(missing_names)
            )
        requested_set = set(requested_canonicals)
        candidate_rows = {
            canonical: candidate_rows.get(canonical, [])
            for canonical in requested_canonicals
            if canonical in candidate_rows or canonical in compound_templates
        }

    if not candidate_rows:
        print("No candidate compound packages found. No changes applied.")
        return

    client = None if args.no_npclassifier else NPClassifierClient()
    approved_logs: List[Dict[str, Any]] = []
    skipped_packages: List[str] = []
    apply_all_remaining = bool(args.yes)

    ordered_candidates = sorted(candidate_rows.keys())
    for canonical in ordered_candidates:
        template = compound_templates[canonical]
        candidate_names = {
            normalize_name(canonical),
            *(normalize_name(name) for name in template.get("alternative_names") or []),
        }
        db_alias_lookup = build_db_compound_alias_lookup(db["compounds"])
        compound_change_key, compound_before, _compound_after = compound_before_after(db, canonical, template)
        replacement_names = set(candidate_names)
        if compound_change_key:
            replacement_names.add(normalize_name(compound_change_key))

        reaction_changes: List[Dict[str, Any]] = []
        unresolved_tokens: List[str] = []
        for activity_row in candidate_rows.get(canonical, []):
            reaction_change, referenced_canonicals, row_unresolved = build_reaction_change(
                db,
                activity_row,
                canonical,
                candidate_names,
                replacement_names,
            )
            unresolved_tokens.extend(row_unresolved)
            if reaction_change and canonical in referenced_canonicals:
                reaction_changes.append(reaction_change)

        affected_row_indices = [change["row_index"] for change in reaction_changes]
        compound_change_key, compound_before, compound_after = compound_before_after(
            db,
            canonical,
            template,
            extra_row_indices=affected_row_indices,
        )
        compound_diff = diff_fields(compound_before, compound_after, COMPOUND_DIFF_FIELDS)

        preview_db = copy.deepcopy(db)
        compound_change = {
            "canonical": canonical,
            "existing_key": compound_change_key,
            "before": compound_before,
            "after": compound_after,
            "diff": compound_diff,
        }
        if compound_diff:
            apply_compound_change(preview_db, canonical, compound_change)
        for reaction_change in reaction_changes:
            apply_reaction_change(preview_db, reaction_change)

        enzyme_changes: List[Dict[str, Any]] = []
        for enzyme_id in collect_changed_enzymes(reaction_changes):
            enzyme_change = build_enzyme_change(preview_db, enzyme_id, canonical, replacement_names)
            if enzyme_change:
                enzyme_changes.append(enzyme_change)

        if not compound_diff and not reaction_changes and not enzyme_changes:
            continue

        npclassifier_needed = needs_npclassifier_fetch(compound_after)
        print_compound_preview(
            canonical,
            compound_change,
            reaction_changes,
            enzyme_changes,
            npclassifier_needed=npclassifier_needed,
            unresolved_tokens=unresolved_tokens,
        )

        if apply_all_remaining:
            decision = "y"
        else:
            decision = prompt_for_approval(canonical)

        if decision == "q":
            skipped_packages.extend(
                candidate for candidate in ordered_candidates[ordered_candidates.index(canonical):]
            )
            break
        if decision == "n":
            skipped_packages.append(canonical)
            continue
        if decision == "a":
            apply_all_remaining = True

        if compound_diff:
            apply_compound_change(db, canonical, compound_change)

        npclassifier_payload = None
        npclassifier_applied = False
        if npclassifier_needed and client is not None:
            smiles = normalize_smiles(compound_after.get("smiles", NA))
            npclassifier_payload = client.classify(smiles)
            if npclassifier_payload is not None:
                db["compounds"][canonical]["npclassifier"] = sanitize_npclassifier(npclassifier_payload)
                npclassifier_applied = True
            else:
                db["compounds"][canonical]["npclassifier"] = {}

        for reaction_change in reaction_changes:
            apply_reaction_change(db, reaction_change)

        applied_enzyme_logs: List[Dict[str, Any]] = []
        for enzyme_change in enzyme_changes:
            refreshed_change = build_enzyme_change(db, enzyme_change["enzyme_id"], canonical, replacement_names)
            if refreshed_change:
                apply_enzyme_change(db, refreshed_change)
                applied_enzyme_logs.append({
                    "enzyme_id": refreshed_change["enzyme_id"],
                    "diff": refreshed_change["diff"],
                })

        approved_logs.append({
            "compound": canonical,
            "matched_existing_key": compound_change_key,
            "compound_diff": compound_diff,
            "npclassifier_requested": bool(npclassifier_needed and client is not None),
            "npclassifier_applied": npclassifier_applied,
            "npclassifier_result": sanitize_npclassifier(npclassifier_payload),
            "reaction_updates": [
                {
                    "row_index": change["row_index"],
                    "enzyme_id": change["enzyme_id"],
                    "diff": change["diff"],
                }
                for change in reaction_changes
            ],
            "enzyme_updates": applied_enzyme_logs,
            "unresolved_activity_tokens": sorted_unique(unresolved_tokens),
        })

    append_log_entry(db, args, approved_logs, skipped_packages)
    write_json_atomic(args.output_json, db)

    print()
    print(f"Approved packages: {len(approved_logs)}")
    if approved_logs:
        print("Applied compounds:", ", ".join(entry["compound"] for entry in approved_logs))
    print(f"Skipped packages: {len(skipped_packages)}")
    if skipped_packages:
        print("Skipped compounds:", ", ".join(skipped_packages))
    print(f"Updated DB written to {args.output_json}")


if __name__ == "__main__":
    main()
