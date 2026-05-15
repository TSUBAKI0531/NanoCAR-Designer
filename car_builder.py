"""Step 5: CAR construct builder.

Modular assembly of chimeric antigen receptor constructs from
recognition domain, hinge, transmembrane, costimulatory, and
signaling domains.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

DATA_DIR = Path(__file__).parent / "data"


@dataclass
class CARDomain:
    """A single domain component of a CAR construct."""

    name: str
    domain_type: str  # recognition, hinge, TM, costim, signaling
    sequence: str
    source: str = ""
    notes: str = ""
    length: int = 0

    def __post_init__(self) -> None:
        self.length = len(self.sequence)


@dataclass
class CARConstruct:
    """Complete CAR construct with all domains assembled."""

    tag_recognition_domain: CARDomain
    hinge: CARDomain
    transmembrane: CARDomain
    costimulatory: CARDomain
    signaling: CARDomain
    full_sequence: str = ""
    generation: str = "2nd"
    total_length: int = 0
    domain_map: list[dict] = field(default_factory=list)

    def __post_init__(self) -> None:
        self.full_sequence = (
            self.tag_recognition_domain.sequence
            + self.hinge.sequence
            + self.transmembrane.sequence
            + self.costimulatory.sequence
            + self.signaling.sequence
        )
        self.total_length = len(self.full_sequence)
        self._build_domain_map()

    def _build_domain_map(self) -> None:
        """Build positional map of domains for visualization."""
        pos = 0
        self.domain_map = []
        for name, domain in [
            ("Recognition", self.tag_recognition_domain),
            ("Hinge", self.hinge),
            ("Transmembrane", self.transmembrane),
            ("Costimulatory", self.costimulatory),
            ("Signaling (CD3ζ)", self.signaling),
        ]:
            end = pos + domain.length
            self.domain_map.append({
                "domain": name,
                "name": domain.name,
                "start": pos,
                "end": end,
                "length": domain.length,
            })
            pos = end


def load_domain_library() -> dict:
    """Load CAR domain reference sequences from JSON.

    Returns:
        Dict with keys: recognition_domains, hinge_domains,
        transmembrane_domains, costimulatory_domains, signaling_domains.
    """
    path = DATA_DIR / "car_domains.json"
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def get_recognition_domain_for_tag(tag_name: str) -> CARDomain:
    """Get the CAR recognition domain matching a tag system.

    Args:
        tag_name: Tag system name (e.g., '5B9', 'PNE', 'ALFA').

    Returns:
        CARDomain for the tag's recognition partner.

    Raises:
        ValueError: If no recognition domain is found for the tag.
    """
    lib = load_domain_library()
    tag_to_key = {
        "5B9": "anti_5B9_scFv",
        "PNE": "anti_GCN4_scFv",
        "FITC": "anti_FITC_scFv",
        "ALFA": "NbALFA",
        "LeucineZipper": "leucine_zipper_BZip",
    }

    key = tag_to_key.get(tag_name)
    if key is None or key not in lib["recognition_domains"]:
        raise ValueError(f"No recognition domain found for tag: {tag_name}")

    entry = lib["recognition_domains"][key]
    return CARDomain(
        name=entry["name"],
        domain_type="recognition",
        sequence=entry["sequence"],
        source=entry.get("source", ""),
        notes=f"Targets: {entry.get('target_tag', '')} | Type: {entry.get('type', '')}",
    )


def get_domain(category: str, domain_key: str) -> CARDomain:
    """Get a specific domain by category and key.

    Args:
        category: Domain category (e.g., 'hinge_domains').
        domain_key: Key within the category (e.g., 'CD8a').

    Returns:
        CARDomain instance.

    Raises:
        ValueError: If the domain is not found.
    """
    lib = load_domain_library()
    cat = lib.get(category, {})
    entry = cat.get(domain_key)
    if entry is None:
        raise ValueError(f"Domain not found: {category}/{domain_key}")

    type_map = {
        "hinge_domains": "hinge",
        "transmembrane_domains": "TM",
        "costimulatory_domains": "costim",
        "signaling_domains": "signaling",
    }

    return CARDomain(
        name=entry["name"],
        domain_type=type_map.get(category, "unknown"),
        sequence=entry["sequence"],
        source=entry.get("source", ""),
        notes=entry.get("notes", ""),
    )


def list_available_domains() -> dict[str, list[str]]:
    """List all available domain options by category.

    Returns:
        Dict mapping category names to lists of domain keys.
    """
    lib = load_domain_library()
    result = {}
    for category in [
        "hinge_domains",
        "transmembrane_domains",
        "costimulatory_domains",
        "signaling_domains",
    ]:
        result[category] = list(lib.get(category, {}).keys())
    return result


def build_car(
    tag_name: str,
    hinge_key: str = "CD8a",
    tm_key: str = "CD8a_TM",
    costim_key: str = "4_1BB",
    signaling_key: str = "CD3z",
) -> CARConstruct:
    """Build a complete CAR construct from modular domains.

    Args:
        tag_name: Tag system name for recognition domain selection.
        hinge_key: Hinge domain key.
        tm_key: Transmembrane domain key.
        costim_key: Costimulatory domain key.
        signaling_key: Signaling domain key.

    Returns:
        CARConstruct with all domains assembled.
    """
    recognition = get_recognition_domain_for_tag(tag_name)
    hinge = get_domain("hinge_domains", hinge_key)
    tm = get_domain("transmembrane_domains", tm_key)
    costim = get_domain("costimulatory_domains", costim_key)
    signaling = get_domain("signaling_domains", signaling_key)

    return CARConstruct(
        tag_recognition_domain=recognition,
        hinge=hinge,
        transmembrane=tm,
        costimulatory=costim,
        signaling=signaling,
    )


def get_car_summary(car: CARConstruct) -> dict:
    """Generate a summary of the CAR construct for display.

    Args:
        car: CARConstruct instance.

    Returns:
        Dict with summary information.
    """
    return {
        "total_length_aa": car.total_length,
        "generation": car.generation,
        "domains": [
            {
                "Domain": d["domain"],
                "Component": d["name"],
                "Length (aa)": d["length"],
                "Position": f"{d['start']+1}-{d['end']}",
            }
            for d in car.domain_map
        ],
        "recognition_type": car.tag_recognition_domain.notes,
        "costimulatory": car.costimulatory.name,
    }
