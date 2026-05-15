# 🧬 NanoCAR Designer

**Modular CAR-T Cell Therapy Design Platform**

Computational design tool for nanobody × tag recognition modular CAR-T constructs. Supports UniCAR, sCAR, SUPRA CAR, Anti-FITC CAR, and ALFA-CAR platforms.

## Pipeline

| Step | Module | Description |
|------|--------|-------------|
| 1 | Biomarker Input | UniProt API lookup or manual sequence entry |
| 2 | Nanobody Design | DB search (RCSB PDB) + template CDR grafting |
| 3 | Tag Selection | Multi-criteria comparison with radar visualization |
| 4 | Fusion Design | Nb-linker-Tag assembly with property calculation |
| 5 | CAR Build | Modular domain selection (hinge/TM/costim/signal) |
| 6 | 3D Structure | ESMFold prediction with py3Dmol visualization |
| 7 | PK Simulation | One-compartment model with ON/OFF dynamics |

## Quick Start

```bash
pip install -r requirements.txt
streamlit run app.py
```

## Architecture

```
NanoCAR-Designer/
├── app.py                    # Streamlit UI (7-step pipeline)
└── src/
    ├── biomarker.py          # UniProt integration
    ├── nanobody.py           # VHH search & CDR grafting
    ├── tag_system.py         # Tag comparison & recommendation
    ├── fusion_designer.py    # Fusion protein design
    ├── car_builder.py        # CAR domain assembly
    ├── structure_viewer.py   # 3D visualization
    ├── pk_simulator.py       # PK/safety simulation
    └── data/                 # Reference databases (JSON)
```

## Tag Systems Supported

- **5B9** (UniCAR) — La/SS-B epitope, Phase I clinical
- **PNE** (sCAR) — GCN4-derived neo-epitope
- **FITC** (Anti-FITC CAR) — Small molecule conjugate
- **ALFA** (ALFA-CAR) — NbALFA-recognized helical peptide
- **Leucine Zipper** (SUPRA CAR) — Tunable coiled-coil system

## References

Based on comprehensive review of nanobody × tag recognition modular CAR-T cell therapy literature (2016-2025).
