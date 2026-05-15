"""Step 4: Nanobody-Tag fusion construct designer.

Assembles nanobody + linker + tag into a complete fusion protein,
calculates physicochemical properties, and generates DNA sequence.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

DATA_DIR = Path(__file__).parent / "data"

# Average amino acid molecular weights (Da)
AA_WEIGHTS: dict[str, float] = {
    "A": 89.1, "R": 174.2, "N": 132.1, "D": 133.1, "C": 121.2,
    "E": 147.1, "Q": 146.2, "G": 75.0, "H": 155.2, "I": 131.2,
    "L": 131.2, "K": 146.2, "M": 149.2, "F": 165.2, "P": 115.1,
    "S": 105.1, "T": 119.1, "W": 204.2, "Y": 181.2, "V": 117.1,
}

# pKa values for pI calculation
PK_VALUES: dict[str, float] = {
    "D": 3.65, "E": 4.25, "C": 8.18, "Y": 10.07,
    "H": 6.00, "K": 10.53, "R": 12.48,
    "N_term": 9.69, "C_term": 2.34,
}

# Codon table for reverse translation (E.coli optimized)
CODON_TABLE_ECOLI: dict[str, str] = {
    "A": "GCG", "R": "CGT", "N": "AAC", "D": "GAT", "C": "TGC",
    "E": "GAA", "Q": "CAG", "G": "GGT", "H": "CAT", "I": "ATT",
    "L": "CTG", "K": "AAA", "M": "ATG", "F": "TTC", "P": "CCG",
    "S": "AGC", "T": "ACC", "W": "TGG", "Y": "TAT", "V": "GTG",
    "*": "TAA",
}

CODON_TABLE_HUMAN: dict[str, str] = {
    "A": "GCC", "R": "CGG", "N": "AAC", "D": "GAC", "C": "TGC",
    "E": "GAG", "Q": "CAG", "G": "GGC", "H": "CAC", "I": "ATC",
    "L": "CTG", "K": "AAG", "M": "ATG", "F": "TTC", "P": "CCC",
    "S": "AGC", "T": "ACC", "W": "TGG", "Y": "TAC", "V": "GTG",
    "*": "TGA",
}


@dataclass
class FusionConstruct:
    """Container for a Nb-Tag fusion protein construct."""

    nanobody_seq: str
    linker_seq: str
    tag_seq: str
    full_sequence: str
    orientation: str  # "Nb-linker-Tag" | "Tag-linker-Nb"
    molecular_weight: float = 0.0
    isoelectric_point: float = 0.0
    instability_index: float = 0.0
    gravy: float = 0.0
    dna_sequence: str = ""
    total_length: int = 0

    def __post_init__(self) -> None:
        self.total_length = len(self.full_sequence)


def load_linker_library() -> list[dict]:
    """Load available linker sequences from JSON data file.

    Returns:
        List of linker dicts with sequence and properties.
    """
    path = DATA_DIR / "linker_library.json"
    with open(path, encoding="utf-8") as f:
        data = json.load(f)
    return data["linkers"]


def design_fusion(
    nanobody_seq: str,
    tag_seq: str,
    linker_id: str,
    orientation: str = "Nb-linker-Tag",
    codon_organism: str = "human",
) -> FusionConstruct:
    """Design a Nb-Tag fusion construct.

    Args:
        nanobody_seq: Nanobody amino acid sequence.
        tag_seq: Tag amino acid sequence.
        linker_id: Linker ID from the linker library.
        orientation: 'Nb-linker-Tag' or 'Tag-linker-Nb'.
        codon_organism: 'human' or 'ecoli' for reverse translation.

    Returns:
        FusionConstruct with full sequence and properties.

    Raises:
        ValueError: If linker_id is not found.
    """
    linkers = load_linker_library()
    linker = None
    for ln in linkers:
        if ln["id"] == linker_id:
            linker = ln
            break
    if linker is None:
        raise ValueError(f"Linker not found: {linker_id}")

    linker_seq = linker["sequence"]

    if orientation == "Nb-linker-Tag":
        full_seq = nanobody_seq + linker_seq + tag_seq
    else:
        full_seq = tag_seq + linker_seq + nanobody_seq

    props = calculate_properties(full_seq)
    dna = reverse_translate(full_seq, codon_organism)

    return FusionConstruct(
        nanobody_seq=nanobody_seq,
        linker_seq=linker_seq,
        tag_seq=tag_seq,
        full_sequence=full_seq,
        orientation=orientation,
        molecular_weight=props["molecular_weight"],
        isoelectric_point=props["isoelectric_point"],
        instability_index=props["instability_index"],
        gravy=props["gravy"],
        dna_sequence=dna,
    )


def calculate_properties(sequence: str) -> dict:
    """Calculate physicochemical properties of a protein sequence.

    Args:
        sequence: Amino acid sequence.

    Returns:
        Dict with molecular_weight, isoelectric_point,
        instability_index, and gravy.
    """
    sequence = sequence.upper()
    mw = _calculate_mw(sequence)
    pi = _calculate_pi(sequence)
    ii = _calculate_instability_index(sequence)
    gravy = _calculate_gravy(sequence)

    return {
        "molecular_weight": round(mw, 1),
        "isoelectric_point": round(pi, 2),
        "instability_index": round(ii, 2),
        "gravy": round(gravy, 3),
    }


def _calculate_mw(sequence: str) -> float:
    """Calculate molecular weight in Daltons."""
    water = 18.02
    mw = sum(AA_WEIGHTS.get(aa, 0) for aa in sequence)
    mw -= water * (len(sequence) - 1)
    return mw


def _calculate_pi(sequence: str) -> float:
    """Calculate isoelectric point using bisection method."""
    pos_aa = {"K": 0, "R": 0, "H": 0}
    neg_aa = {"D": 0, "E": 0, "C": 0, "Y": 0}
    for aa in sequence:
        if aa in pos_aa:
            pos_aa[aa] += 1
        if aa in neg_aa:
            neg_aa[aa] += 1

    def _charge_at_ph(ph: float) -> float:
        pos_charge = 10 ** (PK_VALUES["N_term"] - ph) / (
            1 + 10 ** (PK_VALUES["N_term"] - ph)
        )
        for aa, count in pos_aa.items():
            if count > 0:
                pos_charge += count * 10 ** (PK_VALUES[aa] - ph) / (
                    1 + 10 ** (PK_VALUES[aa] - ph)
                )
        neg_charge = 10 ** (ph - PK_VALUES["C_term"]) / (
            1 + 10 ** (ph - PK_VALUES["C_term"])
        )
        for aa, count in neg_aa.items():
            if count > 0:
                neg_charge += count * 10 ** (ph - PK_VALUES[aa]) / (
                    1 + 10 ** (ph - PK_VALUES[aa])
                )
        return pos_charge - neg_charge

    low, high = 0.0, 14.0
    for _ in range(100):
        mid = (low + high) / 2
        if _charge_at_ph(mid) > 0:
            low = mid
        else:
            high = mid
    return (low + high) / 2


def _calculate_instability_index(sequence: str) -> float:
    """Calculate instability index (Guruprasad et al., 1990).

    Values > 40 suggest the protein may be unstable in vivo.
    """
    # Simplified DIWV weight table (subset of key dipeptides)
    diwv: dict[str, float] = {
        "WW": 1.0, "WC": 1.0, "WM": 24.68, "CK": 1.0, "CC": 1.0,
        "EE": 1.0, "EG": 1.0, "EK": 1.0, "DD": 1.0, "DG": 1.0,
        "GG": 1.0, "GE": 1.0, "KK": 1.0, "KE": 1.0,
    }
    if len(sequence) < 2:
        return 0.0
    total = 0.0
    for i in range(len(sequence) - 1):
        dipep = sequence[i : i + 2]
        total += diwv.get(dipep, 1.0)
    return (10.0 / len(sequence)) * total


def _calculate_gravy(sequence: str) -> float:
    """Calculate Grand Average of Hydropathicity (GRAVY).

    Kyte-Doolittle scale. Negative = hydrophilic, Positive = hydrophobic.
    """
    kd_scale: dict[str, float] = {
        "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5,
        "E": -3.5, "Q": -3.5, "G": -0.4, "H": -3.2, "I": 4.5,
        "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8, "P": -1.6,
        "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
    }
    if not sequence:
        return 0.0
    total = sum(kd_scale.get(aa, 0) for aa in sequence)
    return total / len(sequence)


def reverse_translate(protein_seq: str, organism: str = "human") -> str:
    """Reverse-translate protein sequence to codon-optimized DNA.

    Args:
        protein_seq: Amino acid sequence.
        organism: 'human' or 'ecoli' for codon optimization.

    Returns:
        DNA sequence (5' to 3').
    """
    table = CODON_TABLE_HUMAN if organism == "human" else CODON_TABLE_ECOLI
    codons = [table.get(aa, "NNN") for aa in protein_seq.upper()]
    return "".join(codons)
