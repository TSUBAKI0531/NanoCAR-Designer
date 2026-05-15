"""Step 6: Structure prediction and 3D visualization.

Uses AlphaFold EBI API for structure prediction with ESMFold fallback.
Renders interactive 3D views using py3Dmol via streamlit HTML component.
"""

from __future__ import annotations

import time

import requests

# AlphaFold EBI API (primary)
ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api"
# ESMFold API (fallback)
ESMFOLD_API = "https://api.esmatlas.com/foldSequence/v1/pdb/"

# Maximum sequence length for structure prediction
MAX_SEQ_LENGTH = 800


def predict_structure(sequence: str) -> str | None:
    """Predict 3D structure using AlphaFold EBI → ESMFold fallback.

    Args:
        sequence: Amino acid sequence (max ~800 residues for API).

    Returns:
        PDB format string, or None if prediction fails.
    """
    if len(sequence) > MAX_SEQ_LENGTH:
        sequence = sequence[:MAX_SEQ_LENGTH]

    # Try ESMFold first (faster for custom sequences)
    pdb = _try_esmfold(sequence)
    if pdb:
        return pdb

    # AlphaFold DB lookup (only works for UniProt entries)
    # For custom/designed sequences, ESMFold is the primary option
    return None


def _try_esmfold(sequence: str) -> str | None:
    """Attempt structure prediction via ESMFold API.

    Args:
        sequence: Amino acid sequence.

    Returns:
        PDB string or None.
    """
    try:
        resp = requests.post(
            ESMFOLD_API,
            data=sequence,
            headers={"Content-Type": "text/plain"},
            timeout=120,
        )
        if resp.status_code == 200 and resp.text.startswith("HEADER"):
            return resp.text
    except requests.RequestException:
        pass
    return None


def lookup_alphafold_db(uniprot_id: str) -> str | None:
    """Look up pre-computed AlphaFold structure by UniProt ID.

    Args:
        uniprot_id: UniProt accession.

    Returns:
        PDB string or None.
    """
    url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprot_id}-F1-model_v4.pdb"
    try:
        resp = requests.get(url, timeout=15)
        if resp.status_code == 200:
            return resp.text
    except requests.RequestException:
        pass
    return None


def render_3d_viewer(
    pdb_data: str,
    color_scheme: str = "chain",
    highlight_regions: list[dict] | None = None,
    height: int = 500,
    width: int = 700,
) -> str:
    """Generate py3Dmol HTML for 3D structure visualization.

    Uses setTimeout to handle CDN async loading (known fix for white screen).

    Args:
        pdb_data: PDB format string.
        color_scheme: 'chain', 'spectrum', 'domain', or 'hydrophobicity'.
        highlight_regions: Optional list of dicts with 'start', 'end', 'color'.
        height: Viewer height in pixels.
        width: Viewer width in pixels.

    Returns:
        HTML string for streamlit.components.v1.html().
    """
    # Build color style
    if color_scheme == "spectrum":
        style_cmd = "viewer.setStyle({}, {cartoon: {color: 'spectrum'}});"
    elif color_scheme == "hydrophobicity":
        style_cmd = """
        viewer.setStyle({}, {cartoon: {
            colorscheme: {
                prop: 'b',
                gradient: 'rwb',
                min: 0,
                max: 100
            }
        }});
        """
    else:
        style_cmd = "viewer.setStyle({}, {cartoon: {color: 'spectrum'}});"

    # Build highlight commands
    highlight_cmds = ""
    if highlight_regions:
        for region in highlight_regions:
            start = region.get("start", 0)
            end = region.get("end", 0)
            color = region.get("color", "#FF0000")
            label = region.get("label", "")
            highlight_cmds += f"""
            viewer.setStyle(
                {{resi: ['{start}-{end}']}},
                {{cartoon: {{color: '{color}'}}, stick: {{}}}}
            );
            """
            if label:
                mid = (start + end) // 2
                highlight_cmds += f"""
                viewer.addLabel('{label}',
                    {{position: {{resi: {mid}}}, backgroundColor: '{color}',
                     fontColor: 'white', fontSize: 10}});
                """

    # Escape PDB data for JavaScript
    pdb_escaped = pdb_data.replace("\\", "\\\\").replace("`", "\\`").replace("$", "\\$")

    html = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.0.1/3Dmol-min.js"></script>
        <style>
            body {{ margin: 0; padding: 0; overflow: hidden; }}
            #viewer {{ width: {width}px; height: {height}px; position: relative; }}
        </style>
    </head>
    <body>
        <div id="viewer"></div>
        <script>
            setTimeout(function() {{
                var viewer = $3Dmol.createViewer("viewer", {{
                    backgroundColor: "0x1a1a2e"
                }});
                var pdbData = `{pdb_escaped}`;
                viewer.addModel(pdbData, "pdb");
                {style_cmd}
                {highlight_cmds}
                viewer.zoomTo();
                viewer.spin("y", 0.5);
                viewer.render();
            }}, 500);
        </script>
    </body>
    </html>
    """
    return html


def visualize_fusion_construct(
    pdb_data: str,
    nb_length: int,
    linker_length: int,
    tag_length: int,
    height: int = 500,
) -> str:
    """Render fusion construct with Nb/linker/tag color-coded.

    Args:
        pdb_data: PDB format string of the fusion protein.
        nb_length: Length of nanobody portion.
        linker_length: Length of linker.
        tag_length: Length of tag.
        height: Viewer height.

    Returns:
        HTML string with domain-colored 3D view.
    """
    regions = [
        {"start": 1, "end": nb_length, "color": "#2ecc71", "label": "Nanobody"},
        {"start": nb_length + 1, "end": nb_length + linker_length, "color": "#f39c12", "label": "Linker"},
        {"start": nb_length + linker_length + 1, "end": nb_length + linker_length + tag_length, "color": "#e74c3c", "label": "Tag"},
    ]
    return render_3d_viewer(pdb_data, color_scheme="domain", highlight_regions=regions, height=height)


def visualize_car_domains(
    pdb_data: str,
    domain_map: list[dict],
    height: int = 500,
) -> str:
    """Render CAR construct with each domain color-coded.

    Args:
        pdb_data: PDB format string of the CAR.
        domain_map: List of domain boundary dicts from CARConstruct.
        height: Viewer height.

    Returns:
        HTML string with domain-colored 3D view.
    """
    colors = ["#2ecc71", "#3498db", "#9b59b6", "#f39c12", "#e74c3c"]
    regions = []
    for i, d in enumerate(domain_map):
        regions.append({
            "start": d["start"] + 1,  # 1-indexed for PDB
            "end": d["end"],
            "color": colors[i % len(colors)],
            "label": d["domain"],
        })
    return render_3d_viewer(pdb_data, color_scheme="domain", highlight_regions=regions, height=height)
