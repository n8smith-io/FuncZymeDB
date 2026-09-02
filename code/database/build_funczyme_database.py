#!/usr/bin/env python3
"""
build_funczyme_database.py

Builds an enzyme, compound, and sequence database with optional sequence handling.

This script performs the following steps:
  1. Read a standardized compound TSV file (with COMPOUND_TYPE for donor, acceptor, product).
  2. Read an enzyme activity TSV file (with SUBSTRATE, DONOR, ACCEPTOR, PRODUCT columns).
  3. Validate each record's declared KINGDOM against the species lineage table
     (--species_lineages); mismatches abort the build unless --skip_kingdom_check.
  4. Build or update a JSON database of enzymes, reactions, and compounds.
  5. Initialize the output FASTA, load any manually curated FASTA sequences, and track them by enzyme.
  6. Optionally resolve manual sequences via BLASTP (--resolve_manual_accessions) to identify 100% identity matches
     for sequences lacking GENBANK_PROT accessions, confirm organism matches, and update enzyme entries while retaining
     the prior identifier as alt_enzyme_id.
  7. Reuse sequences from a previous build's FASTA (--seq_cache), matched by enzyme_id or GENBANK_PROT
     accession, so that only genuinely new records cost an NCBI round-trip.
  8. Fetch sequences from NCBI using updated GENBANK_PROT accessions (--fetch_seqs) and write them to the FASTA.
  9. Append any remaining manual sequences, ensuring fasta_id fields are populated for unresolved accessions.
 10. Optionally merge duplicate sequences in-place in the final FASTA file using seqkit (--merge_seqs).
 11. Output a final JSON file (--output_file) and a final FASTA file (--output_fasta).

Usage example:
python code/database/build_funczyme_database.py \
  --compound_file data/curated/compound_data.tsv \
  --activity_file data/curated/activity_data.tsv \
  --species_lineages data/curated/species_lineages.tsv \
  --output_db results/database/funczymedb.json \
  --output_fasta results/database/sequences.fa \
  --person AB \
  --fetch_seqs --ncbi_email <email> --seq_db protein \
  --manual_fasta data/curated/manual_sequences.fa \
  --merge_seqs --merge_log data/curated/sequence_merge.log
"""

import argparse
import csv
import json
import sys
import subprocess
import re
import os
import time
import hashlib
from collections import defaultdict
from typing import Any, Dict, Set, Optional, List, Tuple
from urllib.error import HTTPError, URLError
from urllib.parse import quote_plus
from urllib.request import Request, urlopen
from Bio import Entrez
from Bio.Blast import NCBIWWW, NCBIXML
from datetime import datetime

NA = "NA"

# =============================================================================
# Define the data model classes
# =============================================================================

class Compound:
    """
    Represents a chemical compound with a canonical name, alternative names,
    chemical identifiers (INCHIKEY, PubChem, ChEBI), and various flags for
    donor/acceptor roles and structural features (aromatic/aliphatic).
    """
    def __init__(self, canonical_name: str, inchi_key: str = NA,
                 pubchem_id: str = NA, chebi_id: str = NA,
                 smiles: str = NA, smiles_verified: bool = False,
                 aromatic: bool = False, aliphatic: bool = False,
                 compound_type: str = NA,
                 compound_notes: str = NA, verified_by: str = NA,
                 row_index: Optional[Set[str]] = None,
                 alternative_names: Optional[Set[str]] = None,
                 npclassifier: Optional[Dict[str, Any]] = None):
        self.canonical_name = canonical_name
        # Store alternative names as a set (canonical name is always included)
        self.alternative_names = alternative_names if alternative_names else {canonical_name}
        self.inchi_key = inchi_key
        self.pubchem_id = pubchem_id
        self.chebi_id = chebi_id
        self.smiles = normalize_smiles(smiles)
        self._smiles_verified = bool(smiles_verified and self.smiles != NA)
        self.aromatic = aromatic
        self.aliphatic = aliphatic
        # Set role flags based on compound_type (case-insensitive)
        self.is_donor = compound_type.strip().lower() == "donor"
        self.is_acceptor = compound_type.strip().lower() == "acceptor"
        self.is_product = compound_type.strip().lower() == "product"
        self.compound_notes = compound_notes
        self.verified_by = verified_by
        self.row_index = row_index if row_index else set()
        self.npclassifier = sanitize_npclassifier(npclassifier)
        self._npclassifier_attempted = bool(self.npclassifier)

    def add_alternative_names(self, names: Set[str]):
        """Add new alternative names."""
        self.alternative_names.update(names)

    def add_row_index(self, index: str):
        """Add a source INDEX (from the original data) to the set."""
        if index and index != NA:
            self.row_index.add(index)

    def update_role(self, compound_type: str):
        """Update role flags (donor/acceptor) based on compound_type string."""
        stype = compound_type.strip().lower()
        if stype == "donor":
            self.is_donor = True
        elif stype == "acceptor":
            self.is_acceptor = True

    def update_smiles(self, smiles: str, is_verified: bool = False):
        """Update the stored SMILES, preferring verified values when available."""
        normalized = normalize_smiles(smiles)
        if normalized == NA:
            return
        if self.smiles == NA or not self.smiles.strip():
            self.smiles = normalized
            self._smiles_verified = bool(is_verified)
            return
        if is_verified and (not self._smiles_verified or self.smiles != normalized):
            self.smiles = normalized
            self._smiles_verified = True

    def set_npclassifier(self, data: Optional[Dict[str, Any]]):
        """Store NPClassifier annotations for this compound."""
        self.npclassifier = sanitize_npclassifier(data)
        self._npclassifier_attempted = True

    def mark_npclassifier_attempted(self):
        """Record that an NPClassifier lookup was attempted even if it failed."""
        self._npclassifier_attempted = True

    @property
    def npclassifier_attempted(self) -> bool:
        return self._npclassifier_attempted

    def has_smiles(self) -> bool:
        """Return True if a usable SMILES string is available."""
        return self.smiles not in (NA, "")

    def to_dict(self) -> Dict:
        """Return a dictionary representation for JSON export."""
        return {
            "canonical_name": self.canonical_name,
            "alternative_names": list(self.alternative_names),
            "inchi_key": self.inchi_key,
            "pubchem_id": self.pubchem_id,
            "chebi_id": self.chebi_id,
            "smiles": self.smiles,
            "aromatic": self.aromatic,
            "aliphatic": self.aliphatic,
            "is_donor": self.is_donor,
            "is_acceptor": self.is_acceptor,
            "is_product": self.is_product,
            "compound_notes": self.compound_notes,
            "verified_by": self.verified_by,
            "row_index": list(self.row_index),
            "npclassifier": self.npclassifier
        }

    def __hash__(self):
        # Use lower-case canonical name for uniqueness in sets.
        return hash(self.canonical_name.lower())

    def __eq__(self, other):
        if not isinstance(other, Compound):
            return False
        return self.canonical_name.lower() == other.canonical_name.lower()


def normalize_smiles(smiles: Optional[str]) -> str:
    """Return a cleaned SMILES string or NA if unavailable."""
    if smiles is None:
        return NA
    value = str(smiles).strip()
    if not value or value.upper() == NA:
        return NA
    return value


def select_smiles(smiles: Optional[str], verified_smiles: Optional[str]) -> Tuple[str, bool]:
    """Choose between verified and original SMILES, preferring verified values."""
    for candidate, is_verified in ((verified_smiles, True), (smiles, False)):
        normalized = normalize_smiles(candidate)
        if normalized != NA:
            return normalized, is_verified
    return NA, False


def sanitize_npclassifier(data: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    """Normalize NPClassifier payloads for consistent storage."""
    if not data:
        return {}
    sanitized: Dict[str, Any] = {}
    for key in ("class_results", "superclass_results", "pathway_results"):
        values = data.get(key)
        if values is None:
            continue
        if not isinstance(values, list):
            values = [values]
        sanitized[key] = [str(v) for v in values if v is not None]
    if "isglycoside" in data:
        sanitized["isglycoside"] = bool(data.get("isglycoside"))
    return sanitized



class NPClassifierClient:
    """Minimal client for the NPClassifier API with throttling and caching."""
    BASE_URL = "https://npclassifier.gnps2.org/classify"

    def __init__(self, min_interval: float = 0.5, timeout: int = 15, max_retries: int = 3):
        self.min_interval = max(min_interval, 0.0)
        self.timeout = timeout
        self.max_retries = max(1, max_retries)
        self._last_request = 0.0
        self._cache: Dict[str, Optional[Dict[str, Any]]] = {}

    def _respect_throttle(self, attempt: int = 0) -> None:
        """Sleep if necessary to respect the configured throttle."""
        elapsed = time.time() - self._last_request
        wait = max(self.min_interval - elapsed, 0.0) + (self.min_interval * attempt if attempt else 0.0)
        if wait > 0:
            time.sleep(wait)

    def classify(self, smiles: str) -> Optional[Dict[str, Any]]:
        """Query NPClassifier for a SMILES string and return parsed annotations."""
        normalized = normalize_smiles(smiles)
        if normalized == NA:
            return None

        if normalized in self._cache:
            return self._cache[normalized]

        last_error: Optional[Exception] = None
        for attempt in range(self.max_retries):
            self._respect_throttle(attempt)
            try:
                request = Request(
                    f"{self.BASE_URL}?smiles={quote_plus(normalized)}",
                    headers={"Accept": "application/json"}
                )
                with urlopen(request, timeout=self.timeout) as response:
                    payload = json.loads(response.read().decode("utf-8"))
            except HTTPError as exc:
                last_error = exc
                self._last_request = time.time()
                if (exc.code == 429 or 500 <= exc.code < 600) and attempt < self.max_retries - 1:
                    continue
                break
            except (URLError, TimeoutError, json.JSONDecodeError, ValueError) as exc:
                last_error = exc
                self._last_request = time.time()
                if attempt < self.max_retries - 1:
                    continue
                break
            except Exception as exc:
                last_error = exc
                self._last_request = time.time()
                if attempt < self.max_retries - 1:
                    continue
                break
            else:
                self._last_request = time.time()
                sanitized = sanitize_npclassifier(payload)
                self._cache[normalized] = sanitized
                return sanitized

        if last_error:
            print(f"Warning: NPClassifier request failed for {normalized}: {last_error}")
        self._cache[normalized] = None
        return None


def annotate_compounds_with_npclassifier(compound_registry: Dict[str, Compound],
                                         client: Optional[NPClassifierClient] = None) -> None:
    """Annotate compounds that have SMILES strings using NPClassifier."""
    if client is None:
        client = NPClassifierClient()

    attempted = 0
    successes = 0
    for compound in compound_registry.values():
        if not compound.has_smiles() or compound.npclassifier_attempted:
            continue
        attempted += 1
        result = client.classify(compound.smiles)
        if result is not None:
            compound.set_npclassifier(result)
            if result:
                successes += 1
        else:
            compound.mark_npclassifier_attempted()

    if attempted:
        print(f"NPClassifier annotations attempted: {attempted}, successful: {successes}")


def reuse_npclassifier_annotations(
    compound_registry: Dict[str, Compound],
    database_path: str,
) -> int:
    """Reuse prior NPClassifier results by canonical or alternative name."""
    with open(database_path, "r", encoding="utf-8") as handle:
        cached_database = json.load(handle)

    cached_by_name: Dict[str, Dict[str, Any]] = {}
    for cached in (cached_database.get("compounds") or {}).values():
        names = [cached.get("canonical_name"), *(cached.get("alternative_names") or [])]
        for name in names:
            cleaned = sanitize_string(name, default="")
            if cleaned:
                cached_by_name[cleaned.lower()] = cached

    reused = 0
    for compound in compound_registry.values():
        cached = cached_by_name.get(compound.canonical_name.lower())
        if not cached or not cached.get("npclassifier"):
            continue
        compound.set_npclassifier(cached["npclassifier"])
        reused += 1
    return reused


class Enzyme:
    """
    Represents an enzyme with its identifying information. The unique enzyme id is computed as
    <ENZYME_COMMON_NAME>_<GENBANK_PROT> (or using SWISSPROT_ID/UNIPROT_ID/ALT_ID/GENBANK_NUC/NA as fallbacks). The original row
    INDEX values are stored in row_index for traceability.
    """
    def __init__(self, common_name: str = NA, full_name: str = NA, organism: str = NA,
                 citation: Optional[Dict[str, str]] = None,
                 sequence_ids: Optional[Dict[str, str]] = None,
                 additional_metadata: Optional[Dict[str, str]] = None,
                 initial_row_index: str = NA,
                 manual_seq: bool = False):
        self.common_name = common_name
        self.full_name = full_name
        self.organism = organism
        self.citations = []
        if citation and citation.get("DOI", NA) != NA:
            self.citations.append({
                "DOI": citation.get("DOI", NA),
                "TITLE": citation.get("TITLE", NA),
                "PMID": citation.get("PMID", NA),
                "row_index": {initial_row_index}
            })
        self.sequence_ids = sequence_ids if sequence_ids else {}
        self.additional_metadata = additional_metadata if additional_metadata else {}
        self.orthology_data = {
            "orthologous_group": NA,
            "panther_family": NA,
            "panther_subfamily": NA,
            "clade": NA
        }
        self.row_index: Set[str] = {initial_row_index}
        if self.sequence_ids.get("GENBANK_PROT", NA) != NA:
            self.enzyme_id = f"{self.common_name}_{self.sequence_ids.get('GENBANK_PROT')}"
        elif self.sequence_ids.get("SWISSPROT_ID", NA) != NA:
            self.enzyme_id = f"{self.common_name}_{self.sequence_ids.get('SWISSPROT_ID')}"
        elif self.sequence_ids.get("UNIPROT_ID", NA) != NA:
            self.enzyme_id = f"{self.common_name}_{self.sequence_ids.get('UNIPROT_ID')}"
        elif self.sequence_ids.get("ALT_ID", NA) != NA:
            self.enzyme_id = f"{self.common_name}_{self.sequence_ids.get('ALT_ID')}"
        elif self.sequence_ids.get("GENBANK_NUC", NA) != NA:
            self.enzyme_id = f"{self.common_name}_{self.sequence_ids.get('GENBANK_NUC')}"
        else:
            self.enzyme_id = f"{self.common_name}_NA"
        self.reactions: Set["Reaction"] = set()
        self.associated_compounds: Set[Compound] = set()
        self.manual_seq = manual_seq
        self.alt_enzyme_id = NA
        # fasta_id Initially, it is set to NA. If this is not NA, the enzyme has a sequence.
        self.fasta_id = NA

    def add_reaction(self, reaction: "Reaction"):
        """Add a Reaction to this enzyme and update associated compounds."""
        self.reactions.add(reaction)
        self.row_index.add(reaction.row_index)  # Use row_index here instead of reaction.index
        for comp in reaction.get_participant_compounds():
            self.associated_compounds.add(comp)

    def to_dict(self) -> Dict:
        """Return a dictionary representation for JSON export."""
        return {
            "enzyme_id": self.enzyme_id,
            "common_name": self.common_name,
            "full_name": self.full_name,
            "organism": self.organism,
            "citations": [ {**cit, "row_index": list(cit["row_index"])} for cit in self.citations ],
            "sequence_ids": self.sequence_ids,
            "additional_metadata": self.additional_metadata,
            "orthology_data": self.orthology_data,  # New line added
            "row_index": list(self.row_index),
            "reactions": [rxn.row_index for rxn in self.reactions],
            "associated_compounds": [comp.canonical_name for comp in self.associated_compounds],
            "manual_seq": self.manual_seq,
            "alt_enzyme_id": self.alt_enzyme_id,
            "fasta_id": self.fasta_id
        }

class Reaction:
    """
    Represents a reaction catalyzed by an enzyme. The original row INDEX is stored in row_index.
    """
    def __init__(self, row_index: str, enzyme: Enzyme,
                 substrates: Optional[Set[Compound]] = None,
                 donors: Optional[Set[Compound]] = None,
                 acceptors: Optional[Set[Compound]] = None,
                 products: Optional[Set[Compound]] = None,
                 reaction_metadata: Optional[Dict[str, str]] = None,
                 rhea_id: str = NA,
                 reaction_specific_id: str = NA):
        self.row_index = row_index  # Unique row index from the source
        self.enzyme = enzyme
        self.substrates = substrates if substrates else set()
        self.donors = donors if donors else set()
        self.acceptors = acceptors if acceptors else set()
        self.products = products if products else set()
        self.reaction_metadata = reaction_metadata if reaction_metadata else {}
        self.rhea_id = rhea_id
        self.reaction_specific_id = reaction_specific_id

    def get_participant_compounds(self) -> Set[Compound]:
        """Return the union of all compound participants."""
        return self.substrates.union(self.donors, self.acceptors, self.products)

    def to_dict(self) -> Dict:
        return {
            "reaction_index": self.row_index,
            "enzyme_id": self.enzyme.enzyme_id,
            "substrates": [comp.canonical_name for comp in self.substrates],
            "donors": [comp.canonical_name for comp in self.donors],
            "acceptors": [comp.canonical_name for comp in self.acceptors],
            "products": [comp.canonical_name for comp in self.products],
            "reaction_metadata": self.reaction_metadata,
            "rhea_id": self.rhea_id,
            "reaction_specific_id": self.reaction_specific_id
        }

    def __hash__(self):
        return hash(self.row_index)

    def __eq__(self, other):
        if not isinstance(other, Reaction):
            return False
        return self.row_index == other.row_index

# =============================================================================
# Functions for parsing and ingesting data
# =============================================================================

def parse_compound_str(compound_str: str, role: Optional[str] = None,
                         index: Optional[str] = NA) -> Set[Compound]:
    """
    Parse a semicolon-delimited substrate list. A fully annotated item has the
    format ``<name> | CID:<pubchem> | <INCHIKEY> | CHEBI:<chebi>``, but only the
    name is required. This permits macromolecules, polymers, unnamed
    intermediates, and other substrates for which small-molecule identifiers do
    not exist.
    Optionally, a role (e.g. donor, acceptor) can be provided to set corresponding flags.
    The source index is added for traceability.
    Returns a set of Compound objects.
    """
    compounds = set()
    cleaned_str = sanitize_string(compound_str, default="")
    if cleaned_str == "" or cleaned_str.upper() == NA:
        compounds.add(Compound(NA))
        return compounds

    for comp_part in cleaned_str.split(";"):
        raw_parts = comp_part.split("|")
        name = sanitize_string(raw_parts[0], default=NA)
        if name == NA:
            continue
        pubchem_field = sanitize_string(
            raw_parts[1] if len(raw_parts) > 1 else NA,
            default=NA,
        )
        if pubchem_field != NA and pubchem_field.upper().startswith("CID:"):
            pubchem_field = pubchem_field.split("CID:", 1)[1]
        pubchem = sanitize_string(pubchem_field, default=NA)
        inchi_key = sanitize_string(
            raw_parts[2] if len(raw_parts) > 2 else NA,
            default=NA,
        )
        chebi_field = sanitize_string(
            raw_parts[3] if len(raw_parts) > 3 else NA,
            default=NA,
        )
        if chebi_field != NA and chebi_field.upper().startswith("CHEBI:"):
            chebi_field = chebi_field.split("CHEBI:", 1)[1]
        chebi = sanitize_string(chebi_field, default=NA)
        compound = Compound(canonical_name=name,
                            inchi_key=inchi_key,
                            pubchem_id=pubchem,
                            chebi_id=chebi)
        compound.add_row_index(sanitize_string(index))
        if role:
            compound.update_role(role)
        compounds.add(compound)
    return compounds


def duplicate_activity_indices(records: List[Dict[str, str]]) -> Dict[str, List[int]]:
    """Return missing/duplicated source INDEX values and one-based data rows."""
    locations: Dict[str, List[int]] = defaultdict(list)
    for row_number, record in enumerate(records, start=1):
        locations[get_sanitized_value(record, "INDEX")].append(row_number)
    return {
        index: row_numbers
        for index, row_numbers in locations.items()
        if index == NA or len(row_numbers) > 1
    }

def import_substrate_file(filename: str) -> Tuple[Dict[str, Compound], Dict[str, Compound]]:
    """
    Reads the standardized substrate TSV file and returns:
      - compound_registry: Dictionary of Compound objects keyed by their canonical name.
      - name_lookup: Dictionary mapping lower-case names (canonical and alternatives) to the Compound.
    """
    compound_registry: Dict[str, Compound] = {}
    name_lookup: Dict[str, Compound] = {}

    try:
        with open(filename, newline="", encoding="utf-8-sig", errors="replace") as tsvfile:
            reader = csv.DictReader(tsvfile, delimiter="\t")
            for row in reader:
                raw_name_field = sanitize_string(row.get("COMPOUND_NAME", ""), default="")
                names_list = []
                for name in raw_name_field.split(";"):
                    cleaned_name = sanitize_string(name, default="")
                    if cleaned_name:
                        names_list.append(cleaned_name)
                if not names_list:
                    continue
                canonical = names_list[0]
                alternative_names = set(names_list)
                aromatic = get_sanitized_value(row, "AROMATIC", default="0") == "1"
                aliphatic = get_sanitized_value(row, "ALIPHATIC", default="0") == "1"
                compound_type = get_sanitized_value(row, "COMPOUND_TYPE")
                smiles_value, smiles_verified = select_smiles(
                    row.get("SMILES", NA),
                    row.get("VERIFIED_SMILES", NA)
                )
                cid_field = get_sanitized_value(row, "CID")
                if cid_field != NA and cid_field.upper().startswith("CID:"):
                    cid_field = cid_field.split("CID:", 1)[1]
                pubchem_id = sanitize_string(cid_field, default=NA)
                chebi_field = get_sanitized_value(row, "CHEBI")
                if chebi_field != NA and chebi_field.upper().startswith("CHEBI:"):
                    chebi_field = chebi_field.split("CHEBI:", 1)[1]
                chebi_id = sanitize_string(chebi_field, default=NA)
                row_id = get_sanitized_value(row, "ROW_ID")
                comp = Compound(
                    canonical_name=canonical,
                    alternative_names=alternative_names,
                    inchi_key=get_sanitized_value(row, "INCHIKEY"),
                    pubchem_id=pubchem_id,
                    chebi_id=chebi_id,
                    smiles=smiles_value,
                    smiles_verified=smiles_verified,
                    aromatic=aromatic,
                    aliphatic=aliphatic,
                    compound_type=compound_type,
                    compound_notes=get_sanitized_value(row, "COMPOUND_NOTES"),
                    verified_by=get_sanitized_value(row, "VERIFIED_BY"),
                    row_index={row_id}
                )
                key_lower = canonical.lower()
                if key_lower in name_lookup:
                    existing_comp = name_lookup[key_lower]
                    existing_comp.add_row_index(row_id)
                    existing_comp.update_role(compound_type)
                    existing_comp.update_smiles(smiles_value, smiles_verified)
                    existing_comp.add_alternative_names(set(names_list))
                    for n in names_list:
                        name_lookup[n.lower()] = existing_comp
                    comp = existing_comp
                else:
                    compound_registry[canonical] = comp
                    for n in names_list:
                        name_lookup[n.lower()] = comp
    except Exception as e:
        print(f"Error reading substrate file {filename}: {e}")
        sys.exit(1)
    return compound_registry, name_lookup

# =============================================================================
# Species/kingdom validation
# =============================================================================

CLADE_FIELDS = tuple(f"CLADE{i}" for i in range(1, 9))
CROSS_PATTERN = re.compile(r"\s[xX×]\s")


def normalize_binomial(value: Optional[str]) -> str:
    """Reduce a species string to 'Genus species' for lineage-table lookups."""
    cleaned = sanitize_string(value, default="")
    if not cleaned or cleaned == NA:
        return ""
    cleaned = re.sub(r"\s+", " ", CROSS_PATTERN.sub(" ", cleaned)).strip()
    tokens = cleaned.split(" ")
    if len(tokens) >= 2:
        return f"{tokens[0]} {tokens[1]}"
    return cleaned


def load_species_clades(tsv_path: str) -> Dict[str, Set[str]]:
    """Map each species in the lineage table to its set of lowercase clade names.

    Both the verbatim BINOMIAL_NAME and its normalized 'Genus species' form are
    indexed so that trinomials (subsp./var.) and hybrids resolve either way.
    """
    clade_map: Dict[str, Set[str]] = {}
    try:
        with open(tsv_path, newline="", encoding="utf-8-sig", errors="replace") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                verbatim = sanitize_string(row.get("BINOMIAL_NAME"), default="")
                if not verbatim or verbatim == NA:
                    continue
                clades = {
                    (row.get(field) or "").strip().lower()
                    for field in CLADE_FIELDS
                }
                clades = {clade for clade in clades if clade and clade != "na"}
                for key in {verbatim, normalize_binomial(verbatim)}:
                    if key:
                        clade_map.setdefault(key, set()).update(clades)
    except Exception as e:
        print(f"Error reading species lineages file {tsv_path}: {e}")
        sys.exit(1)
    return clade_map


def validate_species_kingdoms(species_kingdoms: Dict[Tuple[str, str], Set[str]],
                              clade_map: Dict[str, Set[str]]) -> Tuple[List[str], List[str]]:
    """Cross-reference each row's declared KINGDOM against the lineage table.

    Returns (errors, warnings). An error means the activity record claims a
        kingdom the lineage table does not support or names a species with no
        lineage entry at all. The optional taxonomic identifier preceding the
        clade name is retained as provenance but is not interpreted here.
    """
    errors: List[str] = []
    warnings: List[str] = []
    for (species, kingdom), row_indices in sorted(species_kingdoms.items()):
        rows = ", ".join(sorted(row_indices)[:5])
        if len(row_indices) > 5:
            rows += f", ... ({len(row_indices)} rows)"
        if not species or species == NA:
            errors.append(f"missing SPECIES [rows: {rows}]")
            continue
        clades = clade_map.get(species)
        if clades is None:
            clades = clade_map.get(normalize_binomial(species))
        if clades is None:
            errors.append(
                f"'{species}' has no entry in the species lineage table [rows: {rows}]"
            )
            continue
        parts = [part.strip() for part in kingdom.split("|")] if kingdom and kingdom != NA else []
        declared_name = parts[-1].lower() if parts else ""
        if not declared_name:
            errors.append(f"'{species}' has no KINGDOM value [rows: {rows}]")
            continue
        if declared_name not in clades:
            errors.append(
                f"'{species}' declared KINGDOM '{kingdom}' but lineage table lists "
                f"{sorted(clades)} [rows: {rows}]"
            )
            continue
    return errors, warnings


def ingest_enzyme_record(record: Dict[str, str],
                           enzyme_registry: Dict[str, Enzyme],
                           compound_registry: Dict[str, Compound],
                           name_lookup: Dict[str, Compound],
                           nonmatch_logs: Dict[str, Set[str]]):
    row_index = get_sanitized_value(record, "INDEX")
    common_name = get_sanitized_value(record, "ENZYME_COMMON_NAME")
    full_name = get_sanitized_value(record, "ENZYME_FULL_NAME")
    organism = get_sanitized_value(record, "SPECIES")

    citation = {
        "DOI": get_sanitized_value(record, "DOI"),
        "TITLE": get_sanitized_value(record, "TITLE"),
        "PMID": get_sanitized_value(record, "PMID")
    }

    sequence_ids = {
        "GENBANK_PROT": get_sanitized_value(record, "GENBANK_PROT"),
        "GENBANK_NUC": get_sanitized_value(record, "GENBANK_NUC"),
        "UNIPROT_ID": get_sanitized_value(record, "UNIPROT_ID"),
        "SWISSPROT_ID": get_sanitized_value(record, "SWISSPROT_ID"),
        "ALT_ID": get_sanitized_value(record, "ALT_ID")
    }

    additional_metadata = {
        "GENERAL_ENZYME_FAMILY": get_sanitized_value(record, "GENERAL_ENZYME_FAMILY"),
        "KINGDOM": get_sanitized_value(record, "KINGDOM"),
        "PFAM_DOMAIN": get_sanitized_value(record, "PFAM_DOMAIN"),
        "FAMILY_MEMBERSHIP_REFERENCE": get_sanitized_value(
            record, "FAMILY_MEMBERSHIP_REFERENCE"
        ),
        "FAMILY_MEMBERSHIP_CALL": get_sanitized_value(
            record, "FAMILY_MEMBERSHIP_CALL"
        ),
        "EC_NUMBER": get_sanitized_value(record, "EC_NUMBER")
    }

    manual_seq_val = get_sanitized_value(record, "MANUAL_SEQUENCE", default="0").lower()
    manual_seq = manual_seq_val in {"1", "true", "yes"}

    if sequence_ids.get("GENBANK_PROT", NA) != NA:
        unique_id = f"{common_name}_{sequence_ids.get('GENBANK_PROT')}"
    elif sequence_ids.get("SWISSPROT_ID", NA) != NA:
        unique_id = f"{common_name}_{sequence_ids.get('SWISSPROT_ID')}"
    elif sequence_ids.get("UNIPROT_ID", NA) != NA:
        unique_id = f"{common_name}_{sequence_ids.get('UNIPROT_ID')}"
    elif sequence_ids.get("ALT_ID", NA) != NA:
        unique_id = f"{common_name}_{sequence_ids.get('ALT_ID')}"
    elif sequence_ids.get("GENBANK_NUC", NA) != NA:
        unique_id = f"{common_name}_{sequence_ids.get('GENBANK_NUC')}"
    else:
        unique_id = f"{common_name}_NA"

    unique_id = sanitize_string(unique_id)

    if unique_id in enzyme_registry:
        enzyme = enzyme_registry[unique_id]
        enzyme.row_index.add(row_index)
        updated = False
        for cit in enzyme.citations:
            if cit["DOI"] == citation["DOI"] and cit["TITLE"] == citation["TITLE"] and cit["PMID"] == citation["PMID"]:
                cit["row_index"].add(row_index)
                updated = True
                break
        if not updated and citation.get("DOI", NA) != NA:
            enzyme.citations.append({
                "DOI": citation.get("DOI", NA),
                "TITLE": citation.get("TITLE", NA),
                "PMID": citation.get("PMID", NA),
                "row_index": {row_index}
            })
    else:
        enzyme = Enzyme(common_name, full_name, organism,
                        citation, sequence_ids, additional_metadata,
                        initial_row_index=row_index,
                        manual_seq=manual_seq)
        enzyme_registry[unique_id] = enzyme

    rhea_id = get_sanitized_value(record, "RHEA_ID")

    def match_and_update(compounds: Set[Compound], col_name: str) -> Set[Compound]:
        matched = set()
        for comp in compounds:
            key = comp.canonical_name.lower()
            if key in name_lookup:
                substrate_comp = name_lookup[key]
                substrate_comp.add_row_index(row_index)
                if comp.canonical_name not in substrate_comp.alternative_names:
                    substrate_comp.add_alternative_names({comp.canonical_name})
                if comp.is_donor:
                    substrate_comp.is_donor = True
                if comp.is_acceptor:
                    substrate_comp.is_acceptor = True
                name_lookup[comp.canonical_name.lower()] = substrate_comp
                matched.add(substrate_comp)
                # The two lines of code below this comment may result in a compound being identified as both a product and another compound type
                # This may be appropriate, but I haven't yet made a final determination
                if col_name == "PRODUCT":
                    substrate_comp.is_product = True
            else:
                nonmatch_logs[col_name].add(comp.canonical_name)
                matched.add(comp)
                # any brand‐new compound from the PRODUCT column is a product
                if col_name == "PRODUCT":
                    comp.is_product = True
                compound_registry[comp.canonical_name] = comp
                name_lookup[comp.canonical_name.lower()] = comp
        return matched

    substrates = match_and_update(parse_compound_str(record.get("SUBSTRATE", NA), index=row_index), "SUBSTRATE")
    donors = match_and_update(parse_compound_str(record.get("DONOR", NA), role="donor", index=row_index), "DONOR")
    acceptors = match_and_update(parse_compound_str(record.get("ACCEPTOR", NA), role="acceptor", index=row_index), "ACCEPTOR")
    products = match_and_update(parse_compound_str(record.get("PRODUCT", NA), index=row_index), "PRODUCT")

    all_compounds = substrates.union(donors, acceptors, products)
    sorted_names = sorted(comp.canonical_name for comp in all_compounds if comp.canonical_name != NA)
    reaction_specific_id = "_+_".join(sorted_names) if sorted_names else NA
    reaction_metadata = {"row_index": row_index}
    for field in (
        "SOURCE_DATASET",
        "ACQUISITION_METHOD",
        "CURATION_STATUS",
        "EVIDENCE_TYPE",
        "OTHER_COMMENTS",
        "CURATED_BY",
    ):
        reaction_metadata[field] = get_sanitized_value(record, field)

    reaction = Reaction(
        row_index,
        enzyme,
        substrates=substrates,
        donors=donors,
        acceptors=acceptors,
        products=products,
        reaction_metadata=reaction_metadata,
        rhea_id=rhea_id,
        reaction_specific_id=reaction_specific_id
    )
    enzyme.add_reaction(reaction)

def export_data(enzyme_registry: Dict[str, Enzyme],
                compound_registry: Dict[str, Compound],
                filename: str):
    """
    Exports the enzyme, compound, and reaction registries to a JSON file.
    Reactions are aggregated from each enzyme's reaction set and stored as a dictionary
    keyed by the reaction's row_index.
    """
    reactions = {}
    for enzyme in enzyme_registry.values():
        for rxn in enzyme.reactions:
            reactions[rxn.row_index] = rxn.to_dict()

    data = {
        "enzymes": {eid: enzyme.to_dict() for eid, enzyme in enzyme_registry.items()},
        "compounds": {name: comp.to_dict() for name, comp in compound_registry.items()},
        "reactions": reactions
    }
    try:
        with open(filename, "w") as f:
            json.dump(data, f, indent=2)
        print(f"Data successfully exported to {filename}")
    except Exception as e:
        print(f"Error exporting data to {filename}: {e}")

def update_fasta_ids(enzyme_registry: Dict[str, Enzyme], header_mapping: Dict[str, str]) -> None:
    """
    Update the fasta_id attribute for each enzyme based on the merged header mapping.
    This ensures that all enzymes with identical sequences will have the same fasta_id,
    even if that differs from their original enzyme_id.
    """
    for enzyme in enzyme_registry.values():
        if enzyme.fasta_id != NA:
            # Look up the current fasta_id in the mapping; if found, update it.
            enzyme.fasta_id = header_mapping.get(enzyme.fasta_id, enzyme.fasta_id)


def reuse_fasta_ids_from_database(
    enzyme_registry: Dict[str, Enzyme],
    database_path: str,
    fasta_path: str,
) -> int:
    """Restore exact-sequence alias mappings lost from a deduplicated FASTA."""
    with open(database_path, "r", encoding="utf-8") as handle:
        cached_enzymes = (json.load(handle).get("enzymes") or {})
    fasta_headers = set(index_sequence_cache(fasta_path)[0])

    reused = 0
    for enzyme_id, enzyme in enzyme_registry.items():
        cached_fasta_id = sanitize_string(
            (cached_enzymes.get(enzyme_id) or {}).get("fasta_id"),
        )
        if cached_fasta_id == NA or cached_fasta_id not in fasta_headers:
            continue
        if enzyme.fasta_id != cached_fasta_id:
            enzyme.fasta_id = cached_fasta_id
            reused += 1
    return reused

def sanitize_string(value: Optional[str], default: str = NA) -> str:
    """
    Strip surrounding whitespace, drop non-ASCII characters, and normalize blank entries.
    """
    if value is None:
        return default
    cleaned = str(value).strip()
    if not cleaned:
        return default
    cleaned = cleaned.encode("ascii", "ignore").decode("ascii")
    if not cleaned:
        return default
    return cleaned


def get_sanitized_value(row: Dict[str, str], key: str, default: str = NA) -> str:
    """Fetch a value from a TSV row and sanitize whitespace/non-ASCII characters."""
    return sanitize_string(row.get(key, default), default=default)

# =============================================================================
# Functions for sequence fetching
# =============================================================================

def fetch_fasta(accession: str, db: str) -> Optional[str]:
    """Fetch the FASTA record for a given GenBank accession."""
    try:
        handle = Entrez.efetch(db=db, id=accession, rettype="fasta", retmode="text")
        fasta_data = handle.read()
        handle.close()
        return fasta_data
    except Exception as e:
        print(f"Error fetching {accession}: {e}")
        return None

def fetch_sequences_for_enzymes(enzyme_registry: Dict[str, Enzyme], output_fasta: str, email: str, db: str,
                                failed_log: str = "failed_seq_fetches.log") -> None:
    """
    Fetch FASTA sequences for enzymes using only the GENBANK_PROT accession.
    When a sequence is successfully fetched, update the GENBANK_PROT field with the refined accession,
    and recompute enzyme_id and fasta_id accordingly.
    Manual sequences are eligible once a GENBANK_PROT accession has been assigned.
    Enzymes already satisfied from --seq_cache are skipped, so only genuinely
    new records cost a network round-trip. Appends to output_fasta, which main()
    has already truncated.
    """
    Entrez.email = email
    os.makedirs(os.path.dirname(failed_log) or ".", exist_ok=True)
    with open(output_fasta, "a") as outfile, open(failed_log, "w") as logf:
        for enzyme in enzyme_registry.values():
            # Skip anything already satisfied from the sequence cache.
            if enzyme.fasta_id != NA:
                continue
            # Only fetch if a valid GENBANK_PROT accession is available.
            accession = enzyme.sequence_ids.get("GENBANK_PROT", NA)
            if accession == NA:
                continue  # Skip enzymes without a valid GENBANK_PROT accession.
            print(f"Fetching sequence for {enzyme.enzyme_id} (accession: {accession})")
            result = fetch_fasta(accession, db)
            candidate = accession
            # Attempt to clean candidate if necessary:
            if result is None:
                candidate = candidate.replace("�", "")
                if candidate != accession:
                    print(f"Retrying with � removed: {candidate}")
                    result = fetch_fasta(candidate, db)
            if result is None:
                m = re.search(r'\.(\d+)$', candidate)
                if m and m.group(1) != "1":
                    candidate = re.sub(r'\.(\d+)$', ".1", candidate)
                    print(f"Retrying with non-1 integer replaced: {candidate}")
                    result = fetch_fasta(candidate, db)
            if result is None:
                if not re.search(r'\.\d+$', candidate):
                    candidate = candidate + ".1"
                    print(f"Retrying with .1 appended: {candidate}")
                    result = fetch_fasta(candidate, db)
            if result:
                # Update the enzyme's GENBANK_PROT with the refined accession.
                refined_accession = candidate
                enzyme.sequence_ids["GENBANK_PROT"] = refined_accession
                # Recompute enzyme_id and update fasta_id accordingly.
                enzyme.enzyme_id = sanitize_string(f"{enzyme.common_name}_{refined_accession}")
                enzyme.fasta_id = enzyme.enzyme_id
                # Update FASTA header to match the new enzyme_id.
                lines = result.splitlines()
                new_header = f">{enzyme.enzyme_id}"
                lines[0] = new_header
                new_fasta = "\n".join(lines)
                outfile.write(new_fasta)
                if not new_fasta.endswith("\n"):
                    outfile.write("\n")
            else:
                logf.write(f"{accession}\n")
            #time.sleep(0.01)
    print(f"All sequences have been written to {output_fasta}")
    print(f"Failed fetches have been logged to {failed_log}")

def load_manual_sequences(manual_fasta: str) -> Dict[str, str]:
    """Read a FASTA file of manual sequences into a mapping of header -> sequence."""
    sequences: Dict[str, str] = {}
    with open(manual_fasta, "r", encoding="utf-8", errors="replace") as f:
        header = None
        seq_lines: List[str] = []
        for raw_line in f:
            line = raw_line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    sequences[header] = "".join(seq_lines)
                header = line[1:].strip()
                seq_lines = []
            else:
                seq_lines.append(line.strip())
        if header is not None:
            sequences[header] = "".join(seq_lines)
    return sequences


def diagnose_orphan_manual_header(header: str, enzyme_registry: Dict[str, Enzyme]) -> str:
    """Explain why a manual FASTA header no longer matches an enzyme_id.

    An enzyme_id is <ENZYME_COMMON_NAME>_<identifier>, where the identifier is
    the first available of GENBANK_PROT, SWISSPROT_ID, UNIPROT_ID, ALT_ID or
    GENBANK_NUC. Manual FASTA headers are frozen at curation time, so promoting
    a record to a better identifier (typically adding a GENBANK_PROT) silently
    orphans its header. Report the enzyme that shares the header's common name
    so the drift is actionable rather than just a skip.
    """
    candidates = sorted(
        (enzyme.common_name, enzyme.enzyme_id)
        for enzyme in enzyme_registry.values()
        if enzyme.common_name != NA and header.startswith(f"{enzyme.common_name}_")
    )
    if not candidates:
        return "no enzyme in the activity table shares its common name"
    common_name = candidates[0][0]
    return (f"'{common_name}' now resolves to "
            f"{', '.join(enzyme_id for _, enzyme_id in candidates)}")


def map_manual_sequences_to_enzymes(manual_sequences: Dict[str, str],
                                    enzyme_registry: Dict[str, Enzyme]) -> Tuple[Dict[Enzyme, str], List[str]]:
    """Map manual FASTA headers to Enzyme objects, reporting curation drift.

    Returns (mapped, drift). Two kinds of drift are flagged: a header that
    matches no enzyme_id (the record was promoted to a better identifier since
    the sequence was curated), and a header that matches an enzyme the activity
    table marks MANUAL_SEQUENCE = 0 (a real accession has since superseded the
    curated sequence, leaving the FASTA entry redundant). Neither loses data --
    such sequences arrive from --seq_cache or --fetch_seqs instead -- so these
    are warnings, not build failures.
    """
    mapped: Dict[Enzyme, str] = {}
    drift: List[str] = []
    id_lookup = {enzyme.enzyme_id: enzyme for enzyme in enzyme_registry.values()}
    for header, sequence in manual_sequences.items():
        enzyme = id_lookup.get(header)
        if not enzyme:
            drift.append(
                f"'{header}' matches no enzyme_id; "
                f"{diagnose_orphan_manual_header(header, enzyme_registry)}"
            )
            continue
        if not enzyme.manual_seq:
            drift.append(
                f"'{header}' is marked MANUAL_SEQUENCE = 0 in the activity table; "
                "the curated sequence is redundant with the fetched one"
            )
        mapped[enzyme] = sequence
    return mapped, drift


def index_sequence_cache(cache_fasta: str) -> Tuple[Dict[str, str], Dict[str, str]]:
    """Index a previously built FASTA for reuse, avoiding repeat NCBI fetches.

    Returns (by_enzyme_id, by_accession). Cached headers are enzyme_ids of the
    form <common_name>_<accession>, and both halves may contain underscores
    (e.g. Atha_AT3G47170_NP_190301.2), so the accession cannot be split out
    positionally. Instead every underscore-delimited suffix of each header is
    indexed; an enzyme then looks itself up by its own accession. Suffixes that
    are ambiguous across two different sequences are dropped rather than
    guessed at.
    """
    cached = load_manual_sequences(cache_fasta)
    by_enzyme_id: Dict[str, str] = {}
    by_accession: Dict[str, str] = {}
    ambiguous: Set[str] = set()
    for header, sequence in cached.items():
        sequence = re.sub(r"\s+", "", sequence or "")
        if not header or not sequence:
            continue
        by_enzyme_id[header] = sequence
        tokens = header.split("_")
        for i in range(1, len(tokens)):
            suffix = "_".join(tokens[i:])
            if suffix in by_accession and by_accession[suffix] != sequence:
                ambiguous.add(suffix)
            else:
                by_accession[suffix] = sequence
    for suffix in ambiguous:
        by_accession.pop(suffix, None)
    if ambiguous:
        print(f"Note: {len(ambiguous)} ambiguous accession suffix(es) in {cache_fasta} were not indexed.")
    return by_enzyme_id, by_accession


def apply_sequence_cache(enzyme_registry: Dict[str, Enzyme],
                         by_enzyme_id: Dict[str, str],
                         by_accession: Dict[str, str],
                         output_fasta: str,
                         skip: Optional[Set[str]] = None,
                         wrap: int = 60) -> Tuple[int, int]:
    """Write cached sequences for enzymes that have one, returning (hits, misses).

    Matches on enzyme_id first, then falls back to the GENBANK_PROT accession.
    The accession is the stable half of an enzyme_id -- ENZYME_COMMON_NAME is
    hand-edited in the activity table -- so a renamed enzyme still hits.
    Enzymes named in `skip` (those with a curated manual sequence) are left for
    the manual path, which takes precedence over the cache.
    """
    skip = skip or set()
    hits = 0
    misses = 0
    with open(output_fasta, "a") as out:
        for enzyme in enzyme_registry.values():
            if enzyme.fasta_id != NA or enzyme.enzyme_id in skip:
                continue
            sequence = by_enzyme_id.get(enzyme.enzyme_id)
            if not sequence:
                accession = enzyme.sequence_ids.get("GENBANK_PROT", NA)
                if accession != NA:
                    sequence = by_accession.get(accession)
            if not sequence:
                misses += 1
                continue
            enzyme.fasta_id = enzyme.enzyme_id
            out.write(f">{enzyme.enzyme_id}\n")
            for i in range(0, len(sequence), wrap):
                out.write(sequence[i:i+wrap] + "\n")
            hits += 1
    return hits, misses


def extract_species_from_hit(hit_def: str) -> str:
    """Extract the species name from a BLAST hit definition, if present."""
    if not hit_def:
        return ""
    match = re.search(r"\[(.+?)\]$", hit_def)
    if match:
        return match.group(1)
    return ""


def species_matches(reference_species: str, candidate_species: str) -> bool:
    """Return True when both genus and species tokens from reference are found in the candidate string."""
    if not reference_species or reference_species.upper() == NA or not candidate_species:
        return False
    reference_parts = reference_species.split()
    if len(reference_parts) < 2:
        ref_lower = reference_species.lower()
        return ref_lower in candidate_species.lower()
    genus_token = reference_parts[0].lower()
    species_token = reference_parts[1].lower()
    candidate_lower = candidate_species.lower()
    return genus_token in candidate_lower and species_token in candidate_lower


def resolve_accessions_via_blast(enzyme_registry: Dict[str, Enzyme],
                                 manual_enzyme_sequences: Dict[Enzyme, str],
                                 email: str,
                                 blast_db: str = "nr",
                                 hitlist_size: int = 10,
                                 request_delay: float = 0.34) -> None:
    """Resolve GenBank accessions for manual sequences using sequential BLASTP queries."""
    if not manual_enzyme_sequences:
        return

    Entrez.email = email
    registry_keys = {id(enzyme): key for key, enzyme in enzyme_registry.items()}
    resolved_count = 0
    unresolved: List[str] = []

    cache: Dict[str, Optional[str]] = {}
    last_request_time: Optional[float] = None

    for enzyme, raw_sequence in manual_enzyme_sequences.items():
        existing_accession = enzyme.sequence_ids.get("GENBANK_PROT", NA)
        if existing_accession != NA:
            continue

        sequence = re.sub(r"\s+", "", raw_sequence or "")
        if not sequence:
            print(f"Warning: Manual sequence for {enzyme.enzyme_id} is empty; skipping BLAST.")
            unresolved.append(enzyme.enzyme_id)
            continue

        seq_key = hashlib.md5(sequence.encode("utf-8")).hexdigest()
        accession = cache.get(seq_key, None)

        if seq_key not in cache:
            try:
                if request_delay and last_request_time is not None:
                    elapsed = time.monotonic() - last_request_time
                    if elapsed < request_delay:
                        time.sleep(request_delay - elapsed)

                print(f"Running BLASTP for manual enzyme {enzyme.enzyme_id}...")
                handle = NCBIWWW.qblast(
                    "blastp",
                    blast_db,
                    sequence,
                    hitlist_size=hitlist_size,
                    format_type="XML"
                )
                blast_record = NCBIXML.read(handle)
                handle.close()
                last_request_time = time.monotonic()
            except Exception as e:
                print(f"BLAST failed for {enzyme.enzyme_id}: {e}")
                cache[seq_key] = None
                unresolved.append(enzyme.enzyme_id)
                continue

            accession = None
            for alignment in blast_record.alignments:
                species = extract_species_from_hit(getattr(alignment, "hit_def", ""))
                if not species_matches(enzyme.organism, species):
                    continue
                for hsp in alignment.hsps:
                    identities = getattr(hsp, "identities", 0)
                    align_len = getattr(hsp, "align_length", 0)
                    gaps = getattr(hsp, "gaps", 0)
                    if identities == len(sequence) and align_len == len(sequence) and gaps == 0:
                        accession = alignment.accession
                        break
                if accession:
                    break

            cache[seq_key] = accession

        if accession:
            old_id = enzyme.enzyme_id
            if enzyme.alt_enzyme_id == NA:
                enzyme.alt_enzyme_id = old_id
            enzyme.sequence_ids["GENBANK_PROT"] = accession
            enzyme.enzyme_id = sanitize_string(f"{enzyme.common_name}_{accession}")
            enzyme.fasta_id = NA

            old_key = registry_keys.get(id(enzyme), old_id)
            if old_key != enzyme.enzyme_id:
                enzyme_registry[enzyme.enzyme_id] = enzyme
                if old_key in enzyme_registry:
                    del enzyme_registry[old_key]
                registry_keys[id(enzyme)] = enzyme.enzyme_id

            resolved_count += 1
            print(f"Resolved manual enzyme {old_id} to accession {accession} ({enzyme.enzyme_id}).")
        else:
            unresolved.append(enzyme.enzyme_id)

    if resolved_count:
        print(f"Resolved {resolved_count} manual sequence(s) via BLASTP.")
    if unresolved:
        print(f"Manual sequences without matched accessions: {', '.join(sorted(set(unresolved)))}")


def append_manual_sequences(manual_enzyme_sequences: Dict[Enzyme, str],
                            output_fasta: str,
                            wrap: int = 60) -> None:
    """Append manual sequences that still lack a fasta_id to the output FASTA file."""
    if not manual_enzyme_sequences:
        return

    written = 0
    with open(output_fasta, "a") as out:
        for enzyme, raw_sequence in manual_enzyme_sequences.items():
            if enzyme.fasta_id != NA:
                continue
            sequence = re.sub(r"\s+", "", raw_sequence or "")
            if not sequence:
                print(f"Warning: Manual sequence for {enzyme.enzyme_id} is empty; skipping FASTA append.")
                continue
            header = enzyme.enzyme_id
            enzyme.fasta_id = header
            out.write(f">{header}\n")
            for i in range(0, len(sequence), wrap):
                out.write(sequence[i:i+wrap] + "\n")
            written += 1

    if written:
        print(f"Appended {written} manual sequences to {output_fasta}")

# =============================================================================
# Functions for duplicate sequence merging
# =============================================================================

def get_fasta_tab(input_fasta: str) -> List[str]:
    """
    Use seqkit to convert the FASTA file to a tab-delimited format.
    The output is expected to have two columns: header and sequence.
    """
    try:
        result = subprocess.check_output(["seqkit", "fx2tab", input_fasta],
                                         universal_newlines=True)
        return result.strip().splitlines()
    except subprocess.CalledProcessError as e:
        sys.exit("Error running seqkit: " + str(e))
    except FileNotFoundError:
        sys.exit("seqkit not found. Please install seqkit and ensure it is in your PATH.")

def sort_headers(headers: List[str]) -> List[str]:
    """
    Sort a list of FASTA headers by fewest underscores,
    then by length, then lexicographically.
    """
    # 1) find minimal underscore count
    min_unders = min(h.count("_") for h in headers)
    candidates = [h for h in headers if h.count("_") == min_unders]
    # 2) sort by (length, lexicographic)
    return sorted(candidates, key=lambda x: (len(x), x))

def choose_best_header(headers: List[str], manual_merge_dict: Dict[str, str] = None) -> str:
    """
    Given a list of headers (from the same species group for an identical sequence), choose the best header.
    The process:
      1. Filter candidates by the fewest underscore characters.
      2. Reorder candidates by length (shortest first) and lexicographical order.
      3. Compute a candidate key as "|".join(candidates_sorted).
      4. If a manual merge dictionary is provided and the key exists, return the stored selection.
      5. Otherwise, prompt the user for input.
         If the user input is invalid, default to the first (shortest) candidate.
    Returns the chosen header as a string.
    """
    candidates_sorted = sort_headers(headers)
    candidate_key = "|".join(candidates_sorted)

    if manual_merge_dict is not None and candidate_key in manual_merge_dict:
        chosen = manual_merge_dict[candidate_key]
        print(f"Automatically applying saved merge choice for headers [{candidate_key}]: {chosen}")
        return chosen

    if len(candidates_sorted) == 1:
        return candidates_sorted[0]
    else:
        print("Multiple tied headers found:")
        for idx, cand in enumerate(candidates_sorted, start=1):
            print(f"{idx}: {cand}")
        selection = input("Enter the number of the best header: ")
        try:
            selection_idx = int(selection) - 1
            if selection_idx < 0 or selection_idx >= len(candidates_sorted):
                raise ValueError("Selection index out of range")
            chosen = candidates_sorted[selection_idx]
            manual_log = f"Manual selection for tied headers [{candidate_key}]: {', '.join(candidates_sorted)}. Chosen: {chosen}"
            print(manual_log)
            if manual_merge_dict is not None:
                manual_merge_dict[candidate_key] = chosen
            return chosen
        except Exception as e:
            print("Invalid selection. Defaulting to the shortest candidate.")
            default_choice = candidates_sorted[0]
            default_log = f"Invalid manual selection for [{candidate_key}]; defaulted to: {default_choice}"
            if manual_merge_dict is not None:
                manual_merge_dict[candidate_key] = default_choice
            return default_choice

def write_output(merged: Dict, output_fasta: str, wrap: int = 10000) -> None:
    """
    Write the merged unique sequences into the output FASTA file.
    """
    with open(output_fasta, "w") as f:
        for (seq, species_code), header in merged.items():
            f.write(f">{header}\n")
            for i in range(0, len(seq), wrap):
                f.write(seq[i:i+wrap] + "\n")

def write_log(log_entries: List[str], log_file: str) -> None:
    """
    Write the duplicate merge log entries into a log file.
    """
    with open(log_file, "w") as f:
        for entry in log_entries:
            f.write(entry + "\n\n")

def parse_merge_log(log_file: str) -> Dict[str, str]:
    """
    Parse a merge log file (human-readable) and return a manual_merge_dict
    mapping candidate_key ("|"-joined sorted headers) to chosen header.
    """
    manual = {}
    with open(log_file, "r") as f:
        entry_lines = []
        for line in f:
            line = line.rstrip("\n")
            if line == "":
                # process block
                choices = None
                chosen = None
                for l in entry_lines:
                    if l.strip().startswith("Choices:"):
                        choices = [h.strip() for h in l.split("Choices:")[1].split(",")]
                    if l.strip().startswith("Chosen:"):
                        chosen = l.split("Chosen:")[1].strip()
                if choices and chosen:
                    # Sort choices to match choose_best_header
                    sorted_choices = sort_headers(choices)
                    key = "|".join(sorted_choices)
                    manual[key] = chosen
                entry_lines = []
            else:
                entry_lines.append(line)
        # process last block if no trailing blank line
        if entry_lines:
            choices = None
            chosen = None
            for l in entry_lines:
                if l.strip().startswith("Choices:"):
                    choices = [h.strip() for h in l.split("Choices:")[1].split(",")]
                if l.strip().startswith("Chosen:"):
                    chosen = l.split("Chosen:")[1].strip()
            if choices and chosen:
                # Sort choices to match choose_best_header
                sorted_choices = sort_headers(choices)
                key = "|".join(sorted_choices)
                manual[key] = chosen
    return manual

def process_lines(lines: List[str], manual_merge_dict: Dict[str, str] = None):
    header_mapping = {}  # mapping from original header to chosen header
    seq_dict = defaultdict(lambda: defaultdict(list))

    for line in lines:
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        header, sequence = parts[0], parts[1]
        species_code = header[:4]  # assumes species code is first 4 characters
        seq_dict[sequence][species_code].append(header)

    merged = {}
    log_entries = []
    header_mapping = {}
    for seq, species_groups in seq_dict.items():
        for species_code, headers in species_groups.items():
            if len(headers) > 1:
                chosen = choose_best_header(headers, manual_merge_dict)
                # Build a log entry with lines for sorted header choices and the chosen header
                sorted_choices = sort_headers(headers)
                log_entry = (
                    f"Sequence: {seq}\n"
                    f"  Species group ({species_code}) duplicates:\n"
                    f"    Choices: {', '.join(sorted_choices)}\n"
                    f"    Chosen: {chosen}"
                )
                log_entries.append(log_entry)
                for h in headers:
                    header_mapping[h] = chosen
            else:
                chosen = headers[0]
                header_mapping[chosen] = chosen
            merged[(seq, species_code)] = chosen
    return merged, log_entries, header_mapping

def merge_duplicate_seqs(input_fasta: str, output_fasta: str, log_file: str, manual_merge_dict: Dict[str, str] = None) -> Dict[str, str]:
    lines = get_fasta_tab(input_fasta)
    merged, log_entries, header_mapping = process_lines(lines, manual_merge_dict)
    write_output(merged, output_fasta)
    write_log(log_entries, log_file)
    print(f"Merged FASTA file written to: {output_fasta}")
    print(f"Duplicate merge log written to: {log_file}")
    return header_mapping

# =============================================================================
# Main function to parse arguments and process the files
# =============================================================================

def main():
    parser = argparse.ArgumentParser(description="Build and update an enzyme and compound JSON database, with optional GenBank/manual sequence fetching and duplicate merging.")
    #
    parser.add_argument("--compound_file", required=True,
                        help="Path to the standardized compound TSV input file.")
    parser.add_argument("--activity_file", required=True,
                        help="Path to the enzyme activity TSV input file.")
    parser.add_argument("--output_db", type=str, default="enz_db.json",
                        help="Output JSON database file (default: enz_db.json).")
    parser.add_argument("--output_fasta", type=str, default="sequences.fasta",
                        help="Final output FASTA file for fetched, manual, and/or merged sequences.")
    parser.add_argument("--person", required=True,
                        help="The person who is running this code. Initials or name.")
    # Arguments for fetching sequences from GenBank
    parser.add_argument("--fetch_seqs", action="store_true",
                        help="Fetch sequences using GenBank accessions from enzyme data (GENBANK_PROT field).")
    parser.add_argument("--ncbi_email", type=str, default="",
                        help="Email for NCBI Entrez (required if --fetch_seqs is used).")
    parser.add_argument("--seq_db", type=str, default="protein",
                        help="NCBI database to query (default: protein) for sequence fetching.")
    parser.add_argument("--seq_cache", type=str, default="",
                        help="FASTA from a previous build; sequences are reused by enzyme_id or "
                             "GENBANK_PROT accession so --fetch_seqs only retrieves new records.")
    # Arguments for integrating manually gathered sequences
    parser.add_argument("--manual_fasta", type=str, default="",
                        help="Path to manually constructed FASTA file for enzymes with manual_seq True.")
    parser.add_argument("--resolve_manual_accessions", action="store_true",
                        help="Resolve manual sequences to accessions via BLAST. Default: off (append manual sequences as-is).")
    parser.add_argument("--blast_db", type=str, default="nr",
                        help="NCBI BLAST database to use when resolving manual sequences (default: nr).")
    parser.add_argument("--blast_hitlist_size", type=int, default=10,
                        help="Number of BLAST hits to evaluate when searching for identical proteins (default: 10).")
    parser.add_argument("--blast_request_delay", type=float, default=0.34,
                        help="Delay in seconds between BLAST requests to respect NCBI usage guidelines.")
    parser.add_argument("--blast_workers", type=int, default=3,
                        help="Maximum number of concurrent BLAST queries when resolving manual sequences (default: 3).")
    # Arguments for merging duplicate sequences and logging naming convention choices
    parser.add_argument("--merge_seqs", action="store_true",
                        help="Merge duplicate sequences in the output FASTA file using seqkit.")
    parser.add_argument("--merge_log", type=str, default="merge.log",
                        help="Log file for duplicate merge process.")
    parser.add_argument("--species_lineages", type=str,
                        default="data/curated/species_lineages.tsv",
                        help="Species lineage table used to validate each record's declared KINGDOM.")
    parser.add_argument("--skip_kingdom_check", action="store_true",
                        help="Skip validation of SPECIES against the lineage table's clade columns.")
    parser.add_argument(
        "--npclassifier_cache_db",
        default="",
        help=(
            "Optional prior FuncZymeDB JSON whose NPClassifier annotations are "
            "reused before any network lookup."
        ),
    )
    parser.add_argument(
        "--enzyme_cache_db",
        default="",
        help=(
            "Optional prior FuncZymeDB JSON used to restore fasta_id mappings "
            "for aliases removed from the deduplicated sequence cache."
        ),
    )
    args = parser.parse_args()

    # --- Load previous merge choices from the human-readable merge log ---
    manual_merge_dict = {}
    if args.merge_seqs and os.path.exists(args.merge_log):
        manual_merge_dict = parse_merge_log(args.merge_log)
        print(f"Loaded merge choices from {args.merge_log}")

    # --- Build the database ---
    enzyme_registry: Dict[str, Enzyme] = {}
    compound_registry: Dict[str, Compound] = {}
    name_lookup: Dict[str, Compound] = {}
    nonmatch_logs: Dict[str, Set[str]] = {
        "SUBSTRATE": set(),
        "DONOR": set(),
        "ACCEPTOR": set(),
        "PRODUCT": set()
    }

    print("Importing compound file...")
    comp_reg, name_lookup = import_substrate_file(args.compound_file)
    compound_registry.update(comp_reg)
    print(f"Imported {len(compound_registry)} substrate compounds.")

    print("Processing enzyme activity file...")
    species_kingdoms: Dict[Tuple[str, str], Set[str]] = defaultdict(set)
    try:
        with open(args.activity_file, newline="", encoding="utf-8-sig", errors="replace") as tsvfile:
            records = list(csv.DictReader(tsvfile, delimiter="\t"))
    except Exception as e:
        print(f"Error reading enzyme activity file {args.activity_file}: {e}")
        sys.exit(1)

    duplicate_indices = duplicate_activity_indices(records)
    if duplicate_indices:
        details = "; ".join(
            f"{index}: rows {','.join(map(str, row_numbers))}"
            for index, row_numbers in sorted(duplicate_indices.items())
        )
        print(
            "Error: activity INDEX values must be present and unique; " + details,
            file=sys.stderr,
        )
        sys.exit(1)

    for record in records:
        species_kingdoms[
            (get_sanitized_value(record, "SPECIES"), get_sanitized_value(record, "KINGDOM"))
        ].add(get_sanitized_value(record, "INDEX"))
        ingest_enzyme_record(record, enzyme_registry, compound_registry, name_lookup, nonmatch_logs)
    print(f"Processed {len(enzyme_registry)} enzyme records.")

    # --- Validate declared kingdoms against the species lineage table ---
    if args.skip_kingdom_check:
        print("Skipping species/kingdom validation (--skip_kingdom_check).")
    elif not os.path.exists(args.species_lineages):
        print(f"Warning: species lineages file {args.species_lineages} not found; "
              "skipping species/kingdom validation.")
    else:
        clade_map = load_species_clades(args.species_lineages)
        errors, warnings = validate_species_kingdoms(species_kingdoms, clade_map)
        for warning in warnings:
            print(f"Warning: kingdom taxid mismatch: {warning}")
        if errors:
            print(f"Error: {len(errors)} species/kingdom mismatch(es) against "
                  f"{args.species_lineages}:")
            for error in errors:
                print(f" - {error}")
            print("Fix the activity file or the lineage table, or rerun with "
                  "--skip_kingdom_check to bypass.")
            sys.exit(1)
        print(f"Validated {len(species_kingdoms)} species/kingdom pairing(s) against "
              f"{args.species_lineages}.")

    for col in ["SUBSTRATE", "DONOR", "ACCEPTOR", "PRODUCT"]:
        if nonmatch_logs[col]:
            print(f"Non-matches for {col} column:")
            for nm in sorted(nonmatch_logs[col]):
                print(f" - {nm}")

    if args.npclassifier_cache_db:
        reused = reuse_npclassifier_annotations(
            compound_registry,
            args.npclassifier_cache_db,
        )
        print(
            f"Reused NPClassifier annotations for {reused} compound(s) from "
            f"{args.npclassifier_cache_db}."
        )

    print("Annotating compounds with NPClassifier...")
    annotate_compounds_with_npclassifier(compound_registry)

    # --- Sequence Processing ---
    # 1. Initialize/clear the final FASTA
    with open(args.output_fasta, "w") as f:
        pass

    # 2. Load manual sequences (if provided)
    manual_enzyme_sequences: Dict[Enzyme, str] = {}
    if args.manual_fasta:
        if os.path.exists(args.manual_fasta):
            manual_sequences = load_manual_sequences(args.manual_fasta)
            manual_enzyme_sequences, manual_drift = map_manual_sequences_to_enzymes(
                manual_sequences, enzyme_registry
            )
            if manual_drift:
                print(f"Warning: {len(manual_drift)} stale entry(ies) in {args.manual_fasta}; "
                      "no sequences are lost, but these are no longer doing anything:")
                for entry in manual_drift:
                    print(f" - {entry}")
        else:
            print(f"Warning: manual FASTA file {args.manual_fasta} not found.")

    # 3. Resolve missing GENBANK_PROT accessions for manual sequences via BLASTP (optional)
    if manual_enzyme_sequences:
        unresolved_needing_blast = [
            enzyme for enzyme in manual_enzyme_sequences
            if enzyme.sequence_ids.get("GENBANK_PROT", NA) == NA
        ]
        if unresolved_needing_blast and args.resolve_manual_accessions:
            if not args.ncbi_email:
                sys.exit("Error: --ncbi_email is required to resolve manual sequences via BLAST.")
            resolve_accessions_via_blast(
                enzyme_registry,
                manual_enzyme_sequences,
                args.ncbi_email,
                blast_db=args.blast_db,
                hitlist_size=args.blast_hitlist_size,
                request_delay=args.blast_request_delay
            )
        elif unresolved_needing_blast:
            print(f"Skipping manual accession resolution for {len(unresolved_needing_blast)} sequence(s) (--resolve_manual_accessions not set).")

    # 4. Reuse sequences from a previous build before hitting the network
    cache_hits = 0
    cache_misses = 0
    if args.seq_cache:
        if os.path.exists(args.seq_cache):
            by_enzyme_id, by_accession = index_sequence_cache(args.seq_cache)
            cache_hits, cache_misses = apply_sequence_cache(
                enzyme_registry,
                by_enzyme_id,
                by_accession,
                args.output_fasta,
                skip={enzyme.enzyme_id for enzyme in manual_enzyme_sequences}
            )
            print(f"Sequence cache {args.seq_cache}: {cache_hits} reused, "
                  f"{cache_misses} not cached (of {len(enzyme_registry)} enzymes).")
            if args.fetch_seqs and not cache_misses:
                print("All uncached enzymes resolved; no NCBI fetches required.")
        else:
            print(f"Warning: sequence cache {args.seq_cache} not found; fetching all sequences.")

    # 5. Fetch sequences via GenBank (if requested)
    if args.fetch_seqs:
        if not args.ncbi_email:
            sys.exit("Error: --ncbi_email is required when using --fetch_seqs")
        fetch_sequences_for_enzymes(
            enzyme_registry,
            args.output_fasta,
            args.ncbi_email,
            args.seq_db,
            failed_log=os.path.join(os.path.dirname(args.output_db) or ".", "failed_seq_fetches.log")
        )

    # 6. Append any remaining manual sequences to the FASTA and assign fasta_id
    if manual_enzyme_sequences:
        append_manual_sequences(manual_enzyme_sequences, args.output_fasta)

    # 7. Merge duplicate sequences in-place (if requested)
    if args.merge_seqs:
        header_mapping = merge_duplicate_seqs(
            args.output_fasta, args.output_fasta, args.merge_log, manual_merge_dict
        )
        update_fasta_ids(enzyme_registry, header_mapping)

    if args.enzyme_cache_db:
        reused = reuse_fasta_ids_from_database(
            enzyme_registry,
            args.enzyme_cache_db,
            args.output_fasta,
        )
        print(
            f"Restored {reused} cached exact-sequence alias mapping(s) from "
            f"{args.enzyme_cache_db}."
        )

    # --- Export updated JSON database ---
    export_data(enzyme_registry, compound_registry, args.output_db)

    # Add process metadata for build step
    try:
        with open(args.output_db, 'r') as jf:
            db = json.load(jf)
        db.setdefault('log', {})['build_database'] = {
            'time': datetime.now().isoformat(),
            'person': args.person,
            'input': [
                ('compound_file', args.compound_file),
                ('activity_file', args.activity_file),
                ('seq_cache', args.seq_cache or NA),
                ('npclassifier_cache_db', args.npclassifier_cache_db or NA),
                ('enzyme_cache_db', args.enzyme_cache_db or NA)
            ],
            'output': [
                ('output_db', args.output_db),
                ('output_fasta', args.output_fasta)
            ],
            'seq_cache_hits': cache_hits,
            'seq_cache_misses': cache_misses
        }
        with open(args.output_db, 'w') as jf:
            json.dump(db, jf, indent=2)
    except Exception as e:
        print(f"Warning: failed to write process metadata: {e}")

if __name__ == "__main__":
    main()
