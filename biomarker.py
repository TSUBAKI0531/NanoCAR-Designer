"""Step 1: Biomarker sequence input and annotation retrieval.

Fetches protein information from UniProt REST API or accepts manual sequence input.
Extracts extracellular domains suitable for nanobody targeting.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

import requests

UNIPROT_API = "https://rest.uniprot.org/uniprotkb"


@dataclass
class BiomarkerInfo:
    """Container for biomarker protein information."""

    uniprot_id: str
    name: str
    sequence: str
    gene_name: str = ""
    organism: str = ""
    subcellular_location: str = ""
    domains: list[dict] = field(default_factory=list)
    disease_associations: list[str] = field(default_factory=list)
    expression_tissues: list[str] = field(default_factory=list)
    length: int = 0
    molecular_weight: float = 0.0

    def __post_init__(self) -> None:
        self.length = len(self.sequence)


def fetch_from_uniprot(uniprot_id: str) -> BiomarkerInfo:
    """Fetch protein information from UniProt REST API.

    Args:
        uniprot_id: UniProt accession ID (e.g., 'P00533' for EGFR).

    Returns:
        BiomarkerInfo with populated fields.

    Raises:
        ValueError: If the UniProt ID is invalid or not found.
        ConnectionError: If the API is unreachable.
    """
    uniprot_id = uniprot_id.strip().upper()
    if not re.match(r"^[A-Z0-9]{6,10}$", uniprot_id):
        raise ValueError(f"Invalid UniProt ID format: {uniprot_id}")

    url = f"{UNIPROT_API}/{uniprot_id}.json"
    try:
        resp = requests.get(url, timeout=15)
    except requests.RequestException as exc:
        raise ConnectionError(f"UniProt API connection failed: {exc}") from exc

    if resp.status_code == 404:
        raise ValueError(f"UniProt ID not found: {uniprot_id}")
    resp.raise_for_status()

    data = resp.json()
    return _parse_uniprot_json(uniprot_id, data)


def _parse_uniprot_json(uniprot_id: str, data: dict) -> BiomarkerInfo:
    """Parse UniProt JSON response into BiomarkerInfo."""
    # Protein name
    protein_desc = data.get("proteinDescription", {})
    rec_name = protein_desc.get("recommendedName", {})
    name = rec_name.get("fullName", {}).get("value", "Unknown")

    # Sequence
    seq_data = data.get("sequence", {})
    sequence = seq_data.get("value", "")
    mol_weight = seq_data.get("molWeight", 0)

    # Gene name
    genes = data.get("genes", [{}])
    gene_name = ""
    if genes:
        gene_name = genes[0].get("geneName", {}).get("value", "")

    # Organism
    organism = data.get("organism", {}).get("scientificName", "")

    # Subcellular location
    subcellular = ""
    for comment in data.get("comments", []):
        if comment.get("commentType") == "SUBCELLULAR LOCATION":
            locs = comment.get("subcellularLocations", [])
            loc_names = [
                loc.get("location", {}).get("value", "") for loc in locs
            ]
            subcellular = "; ".join(filter(None, loc_names))
            break

    # Domains (features)
    domains = _extract_domains(data.get("features", []))

    # Disease associations
    diseases = []
    for comment in data.get("comments", []):
        if comment.get("commentType") == "DISEASE":
            disease = comment.get("disease", {})
            disease_name = disease.get("diseaseId", "")
            if disease_name:
                diseases.append(disease_name)

    # Tissue expression
    tissues = []
    for comment in data.get("comments", []):
        if comment.get("commentType") == "TISSUE SPECIFICITY":
            for text_obj in comment.get("texts", []):
                tissues.append(text_obj.get("value", ""))

    return BiomarkerInfo(
        uniprot_id=uniprot_id,
        name=name,
        sequence=sequence,
        gene_name=gene_name,
        organism=organism,
        subcellular_location=subcellular,
        domains=domains,
        disease_associations=diseases,
        expression_tissues=tissues,
        molecular_weight=mol_weight,
    )


def _extract_domains(features: list[dict]) -> list[dict]:
    """Extract domain, topological, and region features."""
    domain_types = {
        "Domain", "Region", "Topological domain",
        "Transmembrane", "Signal peptide",
    }
    domains = []
    for feat in features:
        feat_type = feat.get("type", "")
        if feat_type in domain_types:
            loc = feat.get("location", {})
            start = loc.get("start", {}).get("value", 0)
            end = loc.get("end", {}).get("value", 0)
            desc = feat.get("description", "")
            domains.append({
                "type": feat_type,
                "description": desc,
                "start": start,
                "end": end,
            })
    return domains


def parse_manual_input(sequence: str, name: str = "User input") -> BiomarkerInfo:
    """Validate and create BiomarkerInfo from manual sequence input.

    Args:
        sequence: Amino acid sequence (single-letter code).
        name: Optional protein name.

    Returns:
        BiomarkerInfo with basic fields populated.

    Raises:
        ValueError: If the sequence contains invalid characters.
    """
    sequence = re.sub(r"\s+", "", sequence.upper())
    valid_aa = set("ACDEFGHIKLMNPQRSTVWY")
    invalid = set(sequence) - valid_aa
    if invalid:
        raise ValueError(
            f"Invalid amino acid characters: {', '.join(sorted(invalid))}"
        )
    if len(sequence) < 10:
        raise ValueError("Sequence too short (minimum 10 residues)")

    return BiomarkerInfo(
        uniprot_id="MANUAL",
        name=name,
        sequence=sequence,
    )


def get_extracellular_domains(info: BiomarkerInfo) -> list[dict]:
    """Extract extracellular domains suitable for nanobody targeting.

    Args:
        info: BiomarkerInfo with domain annotations.

    Returns:
        List of extracellular domain dicts with sequence excerpts.
    """
    extracellular = []
    for d in info.domains:
        desc_lower = d.get("description", "").lower()
        dtype = d.get("type", "")
        if dtype == "Topological domain" and "extracellular" in desc_lower:
            start = d["start"] - 1  # 0-indexed
            end = d["end"]
            extracellular.append({
                **d,
                "sequence": info.sequence[start:end],
                "length": end - start,
            })
    return extracellular


def calculate_basic_properties(sequence: str) -> dict:
    """Calculate basic protein properties.

    Args:
        sequence: Amino acid sequence.

    Returns:
        Dict with molecular_weight, length, and amino acid composition.
    """
    # Average amino acid molecular weights
    aa_weights = {
        "A": 89.1, "R": 174.2, "N": 132.1, "D": 133.1, "C": 121.2,
        "E": 147.1, "Q": 146.2, "G": 75.0, "H": 155.2, "I": 131.2,
        "L": 131.2, "K": 146.2, "M": 149.2, "F": 165.2, "P": 115.1,
        "S": 105.1, "T": 119.1, "W": 204.2, "Y": 181.2, "V": 117.1,
    }
    mw = sum(aa_weights.get(aa, 0) for aa in sequence) - 18.02 * (len(sequence) - 1)

    composition = {}
    for aa in sequence:
        composition[aa] = composition.get(aa, 0) + 1

    return {
        "molecular_weight": round(mw, 1),
        "length": len(sequence),
        "composition": composition,
    }
