"""Step 2: Nanobody sequence search and template-based design.

Provides two approaches:
  A) Search known nanobody databases (SAbDab) for existing binders.
  B) CDR grafting onto humanized VHH framework templates.
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from pathlib import Path

import requests

DATA_DIR = Path(__file__).parent / "data"
SABDAB_SEARCH_URL = "https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab/search/"


@dataclass
class NanobodyCandidate:
    """Container for a nanobody sequence and its annotations."""

    sequence: str
    source: str  # "database" | "template"
    cdr1: str = ""
    cdr2: str = ""
    cdr3: str = ""
    framework_regions: list[str] = field(default_factory=list)  # FR1-FR4
    target_antigen: str = ""
    origin_db: str | None = None
    pdb_id: str | None = None
    affinity_reported: str | None = None
    framework_id: str | None = None
    validation: dict = field(default_factory=dict)


def load_framework_templates() -> list[dict]:
    """Load VHH framework templates from JSON data file.

    Returns:
        List of framework template dicts.
    """
    path = DATA_DIR / "vhh_frameworks.json"
    with open(path, encoding="utf-8") as f:
        data = json.load(f)
    return data["frameworks"]


def get_cdr_guidelines() -> dict:
    """Load CDR length guidelines.

    Returns:
        Dict with CDR1/CDR2/CDR3 length ranges and notes.
    """
    path = DATA_DIR / "vhh_frameworks.json"
    with open(path, encoding="utf-8") as f:
        data = json.load(f)
    return data["cdr_length_guidelines"]


def search_known_nanobodies(target_name: str) -> list[NanobodyCandidate]:
    """Search for known nanobodies targeting a given antigen.

    Uses SAbDab API and PDB keyword search as fallback.

    Args:
        target_name: Target antigen name (e.g., 'EGFR', 'HER2').

    Returns:
        List of NanobodyCandidate from database hits.
    """
    candidates = []

    # Strategy 1: RCSB PDB search for VHH structures
    pdb_hits = _search_rcsb_vhh(target_name)
    for hit in pdb_hits:
        candidates.append(NanobodyCandidate(
            sequence=hit.get("sequence", ""),
            source="database",
            target_antigen=target_name,
            origin_db="RCSB PDB",
            pdb_id=hit.get("pdb_id"),
            cdr1=hit.get("cdr1", ""),
            cdr2=hit.get("cdr2", ""),
            cdr3=hit.get("cdr3", ""),
        ))

    return candidates


def _search_rcsb_vhh(target_name: str) -> list[dict]:
    """Search RCSB PDB for VHH/nanobody structures against a target.

    Args:
        target_name: Target protein name.

    Returns:
        List of hit dicts with pdb_id and sequence.
    """
    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "full_text",
                    "parameters": {
                        "value": f"nanobody VHH {target_name}"
                    },
                },
            ],
        },
        "return_type": "polymer_entity",
        "request_options": {"results_content_type": ["experimental"], "return_all_hits": False, "paginate": {"start": 0, "rows": 10}},
    }

    url = "https://search.rcsb.org/rcsbsearch/v2/query"
    try:
        resp = requests.post(url, json=query, timeout=15)
        if resp.status_code != 200:
            return []
        data = resp.json()
    except (requests.RequestException, ValueError):
        return []

    hits = []
    for result in data.get("result_set", [])[:5]:
        entity_id = result.get("identifier", "")
        pdb_id = entity_id.split("_")[0] if "_" in entity_id else entity_id
        seq = _fetch_pdb_sequence(entity_id)
        if seq and _is_likely_vhh(seq):
            hits.append({
                "pdb_id": pdb_id,
                "sequence": seq,
                "cdr1": "",
                "cdr2": "",
                "cdr3": "",
            })

    return hits


def _fetch_pdb_sequence(entity_id: str) -> str:
    """Fetch polymer entity sequence from RCSB.

    Args:
        entity_id: RCSB entity identifier (e.g., '7XYZ_1').

    Returns:
        Amino acid sequence string, or empty string on failure.
    """
    url = f"https://data.rcsb.org/rest/v1/core/polymer_entity/{entity_id.replace('_', '/')}"
    try:
        resp = requests.get(url, timeout=10)
        if resp.status_code != 200:
            return ""
        data = resp.json()
        seq = data.get("entity_poly", {}).get("pdbx_seq_one_letter_code_can", "")
        return seq.replace("\n", "")
    except (requests.RequestException, ValueError):
        return ""


def _is_likely_vhh(sequence: str) -> bool:
    """Heuristic check if a sequence is likely a VHH/nanobody.

    Args:
        sequence: Amino acid sequence.

    Returns:
        True if the sequence has VHH-like characteristics.
    """
    length = len(sequence)
    # VHH are typically 110-140 residues
    if length < 100 or length > 160:
        return False
    # Should contain conserved Cys residues
    if sequence.count("C") < 2:
        return False
    return True


def graft_cdrs(
    framework_id: str,
    cdr1: str,
    cdr2: str,
    cdr3: str,
) -> NanobodyCandidate:
    """Graft CDR sequences onto a VHH framework template.

    Args:
        framework_id: ID of the framework template.
        cdr1: CDR1 amino acid sequence.
        cdr2: CDR2 amino acid sequence.
        cdr3: CDR3 amino acid sequence.

    Returns:
        NanobodyCandidate with assembled full sequence.

    Raises:
        ValueError: If framework_id is not found or CDR sequences are invalid.
    """
    templates = load_framework_templates()
    framework = None
    for tmpl in templates:
        if tmpl["id"] == framework_id:
            framework = tmpl
            break

    if framework is None:
        raise ValueError(f"Framework not found: {framework_id}")

    # Validate CDR sequences
    for name, seq in [("CDR1", cdr1), ("CDR2", cdr2), ("CDR3", cdr3)]:
        _validate_cdr(name, seq)

    # Assemble full sequence: FR1-CDR1-FR2-CDR2-FR3-CDR3-FR4
    full_seq = (
        framework["fr1"]
        + cdr1.upper()
        + framework["fr2"]
        + cdr2.upper()
        + framework["fr3"]
        + cdr3.upper()
        + framework["fr4"]
    )

    return NanobodyCandidate(
        sequence=full_seq,
        source="template",
        cdr1=cdr1.upper(),
        cdr2=cdr2.upper(),
        cdr3=cdr3.upper(),
        framework_regions=[
            framework["fr1"],
            framework["fr2"],
            framework["fr3"],
            framework["fr4"],
        ],
        framework_id=framework_id,
        validation=validate_nanobody(full_seq),
    )


def _validate_cdr(name: str, seq: str) -> None:
    """Validate a CDR sequence.

    Args:
        name: CDR name for error messages.
        seq: CDR amino acid sequence.

    Raises:
        ValueError: If the sequence is empty or contains invalid characters.
    """
    if not seq or not seq.strip():
        raise ValueError(f"{name} sequence cannot be empty")
    seq = seq.strip().upper()
    valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
    invalid = set(seq) - valid_aa
    if invalid:
        raise ValueError(
            f"{name} contains invalid characters: {', '.join(sorted(invalid))}"
        )


def annotate_regions(sequence: str) -> dict:
    """Annotate CDR and framework regions using simple heuristics.

    Uses approximate IMGT-like positioning for VHH sequences.
    This is a simplified annotation; for production use, ANARCI or
    IMGT/DomainGapAlign should be used.

    Args:
        sequence: Full VHH amino acid sequence.

    Returns:
        Dict with estimated FR and CDR boundaries.
    """
    seq_len = len(sequence)
    if seq_len < 100 or seq_len > 160:
        return {"error": "Sequence length outside VHH range (100-160 aa)"}

    # Approximate IMGT positions for VHH
    regions = {
        "FR1": {"start": 0, "end": 25, "sequence": sequence[0:25]},
        "CDR1": {"start": 25, "end": 33, "sequence": sequence[25:33]},
        "FR2": {"start": 33, "end": 49, "sequence": sequence[33:49]},
        "CDR2": {"start": 49, "end": 57, "sequence": sequence[49:57]},
        "FR3": {"start": 57, "end": seq_len - 15, "sequence": sequence[57:seq_len - 15]},
        "CDR3": {"start": seq_len - 15, "end": seq_len - 11, "sequence": sequence[seq_len - 15:seq_len - 11]},
        "FR4": {"start": seq_len - 11, "end": seq_len, "sequence": sequence[seq_len - 11:]},
    }
    return regions


def validate_nanobody(sequence: str) -> dict:
    """Validate VHH sequence characteristics.

    Args:
        sequence: Full VHH amino acid sequence.

    Returns:
        Dict with validation results (passed, warnings, errors).
    """
    results = {"passed": [], "warnings": [], "errors": []}
    seq_len = len(sequence)

    # Length check
    if 110 <= seq_len <= 140:
        results["passed"].append(f"Length OK ({seq_len} aa, typical VHH range)")
    elif 100 <= seq_len <= 160:
        results["warnings"].append(f"Length {seq_len} aa is at the edge of VHH range")
    else:
        results["errors"].append(f"Length {seq_len} aa is outside VHH range (100-160)")

    # Conserved Cys residues
    cys_positions = [i for i, aa in enumerate(sequence) if aa == "C"]
    if len(cys_positions) >= 2:
        results["passed"].append(
            f"Conserved Cys found at positions {cys_positions[:2]}"
        )
    else:
        results["errors"].append("Missing conserved Cys residues for disulfide bond")

    # Conserved Trp (usually around position 36 in VHH)
    if "W" in sequence[30:45]:
        results["passed"].append("Conserved Trp found in FR2 region")
    else:
        results["warnings"].append("No Trp found in expected FR2 region (pos 30-45)")

    # VHH hallmark: hydrophilic residues at positions 42, 49, 50, 52
    # (numbering approximate for different frameworks)
    if seq_len > 52:
        hallmark_pos = sequence[42:53]
        hydrophilic = sum(1 for aa in hallmark_pos if aa in "DERHKQN")
        if hydrophilic >= 2:
            results["passed"].append("VHH hallmark hydrophilic residues detected")
        else:
            results["warnings"].append(
                "Low hydrophilic character at VHH hallmark positions"
            )

    # Overall amino acid validity
    valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
    invalid = set(sequence.upper()) - valid_aa
    if invalid:
        results["errors"].append(f"Invalid amino acids: {', '.join(sorted(invalid))}")
    else:
        results["passed"].append("All amino acids valid")

    return results
