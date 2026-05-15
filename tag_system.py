"""Step 3: Tag system selection and comparison.

Loads tag database and provides comparison, recommendation,
and radar chart data for tag system selection.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

DATA_DIR = Path(__file__).parent / "data"


@dataclass
class TagSystem:
    """Container for a tag system and its properties."""

    name: str
    full_name: str
    sequence: str
    molecular_weight_da: float
    immunogenicity: str
    immunogenicity_score: int
    clinical_stage: str
    clinical_stage_score: int
    pk_half_life_min: float
    clearance_route: str
    recognition_domain_type: str
    recognition_domain_name: str
    platform: str
    origin: str
    pros: list[str] = field(default_factory=list)
    cons: list[str] = field(default_factory=list)
    references: list[str] = field(default_factory=list)
    safety_score: int = 0
    efficacy_score: int = 0
    manufacturability_score: int = 0


def load_tag_database() -> list[TagSystem]:
    """Load all tag systems from JSON data file.

    Returns:
        List of TagSystem dataclass instances.
    """
    path = DATA_DIR / "tag_database.json"
    with open(path, encoding="utf-8") as f:
        data = json.load(f)
    return [TagSystem(**entry) for entry in data]


def get_tag_by_name(name: str) -> TagSystem | None:
    """Retrieve a specific tag system by name.

    Args:
        name: Tag system name (e.g., '5B9', 'PNE').

    Returns:
        TagSystem if found, None otherwise.
    """
    for tag in load_tag_database():
        if tag.name == name:
            return tag
    return None


def compare_tags(selected_names: list[str]) -> dict:
    """Generate comparison data for selected tag systems.

    Args:
        selected_names: List of tag names to compare.

    Returns:
        Dict with 'tags' list and 'radar_data' for visualization.
    """
    all_tags = load_tag_database()
    selected = [t for t in all_tags if t.name in selected_names]

    radar_data = []
    for tag in selected:
        radar_data.append({
            "name": tag.name,
            "Safety": tag.safety_score,
            "Efficacy": tag.efficacy_score,
            "Clinical readiness": tag.clinical_stage_score,
            "Low immunogenicity": tag.immunogenicity_score,
            "Manufacturability": tag.manufacturability_score,
        })

    comparison_table = []
    for tag in selected:
        comparison_table.append({
            "Tag": tag.name,
            "Platform": tag.platform,
            "Size (Da)": tag.molecular_weight_da,
            "Immunogenicity": tag.immunogenicity,
            "Clinical stage": tag.clinical_stage,
            "Half-life (min)": tag.pk_half_life_min,
            "CAR recognition": tag.recognition_domain_type,
        })

    return {
        "tags": selected,
        "radar_data": radar_data,
        "comparison_table": comparison_table,
    }


def recommend_tag(priority: str) -> TagSystem:
    """Recommend a tag system based on priority criteria.

    Args:
        priority: One of 'safety', 'efficacy', 'clinical_readiness',
                  'manufacturability'.

    Returns:
        Top-ranked TagSystem for the given priority.
    """
    tags = load_tag_database()

    score_map = {
        "safety": lambda t: t.safety_score + t.immunogenicity_score,
        "efficacy": lambda t: t.efficacy_score + t.clinical_stage_score,
        "clinical_readiness": lambda t: t.clinical_stage_score * 2,
        "manufacturability": lambda t: t.manufacturability_score + t.safety_score,
    }

    key_fn = score_map.get(priority, score_map["safety"])
    ranked = sorted(tags, key=key_fn, reverse=True)
    return ranked[0]


def get_radar_chart_data(tag_names: list[str] | None = None) -> dict:
    """Generate radar chart data for plotly visualization.

    Args:
        tag_names: Optional filter list. If None, all tags are included.

    Returns:
        Dict with 'categories' and 'traces' for plotly radar chart.
    """
    tags = load_tag_database()
    if tag_names:
        tags = [t for t in tags if t.name in tag_names]

    categories = [
        "Safety", "Efficacy", "Clinical\nreadiness",
        "Low\nimmunogenicity", "Manufact-\nurability",
    ]

    traces = []
    for tag in tags:
        values = [
            tag.safety_score,
            tag.efficacy_score,
            tag.clinical_stage_score,
            tag.immunogenicity_score,
            tag.manufacturability_score,
        ]
        # Close the radar polygon
        values.append(values[0])
        cats = categories + [categories[0]]
        traces.append({
            "name": f"{tag.name} ({tag.platform})",
            "categories": cats,
            "values": values,
        })

    return {"categories": categories, "traces": traces}
