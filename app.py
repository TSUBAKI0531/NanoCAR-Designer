"""NanoCAR Designer — Modular CAR-T Design Platform.

Streamlit-based interactive pipeline for designing nanobody × tag
recognition modular CAR-T constructs.

7-step pipeline:
  1. Biomarker sequence input
  2. Nanobody design (DB search + template CDR grafting)
  3. Tag system selection
  4. Nb-Tag fusion construct design
  5. CAR structure design
  6. Structure prediction & 3D visualization
  7. PK / safety simulation
"""

import streamlit as st
import pandas as pd
import plotly.graph_objects as go
import streamlit.components.v1 as components
from dataclasses import asdict

# --- Module imports ---
from biomarker import (
    fetch_from_uniprot,
    parse_manual_input,
    get_extracellular_domains,
    calculate_basic_properties,
)
from nanobody import (
    search_known_nanobodies,
    load_framework_templates,
    get_cdr_guidelines,
    graft_cdrs,
    annotate_regions,
    validate_nanobody,
    NanobodyCandidate,
)
from tag_system import (
    load_tag_database,
    get_tag_by_name,
    compare_tags,
    recommend_tag,
    get_radar_chart_data,
)
from fusion_designer import (
    load_linker_library,
    design_fusion,
    calculate_properties,
)
from car_builder import (
    build_car,
    get_car_summary,
    list_available_domains,
    load_domain_library,
)
from structure_viewer import (
    predict_structure,
    render_3d_viewer,
    visualize_fusion_construct,
    visualize_car_domains,
    lookup_alphafold_db,
)
from pk_simulator import (
    PKParameters,
    simulate_on_off_dynamics,
    generate_pk_report,
)

# ────────────────────────────────────────────
# Page config
# ────────────────────────────────────────────
st.set_page_config(
    page_title="NanoCAR Designer",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ────────────────────────────────────────────
# Custom CSS
# ────────────────────────────────────────────
st.markdown("""
<style>
    @import url('https://fonts.googleapis.com/css2?family=JetBrains+Mono:wght@400;700&family=DM+Sans:wght@400;500;700&display=swap');

    .main .block-container { max-width: 1100px; padding-top: 2rem; }
    h1, h2, h3 { font-family: 'DM Sans', sans-serif; }
    .stCode, code { font-family: 'JetBrains Mono', monospace; font-size: 0.82em; }

    /* Step indicator pills */
    .step-pill {
        display: inline-block;
        padding: 4px 14px;
        border-radius: 20px;
        font-size: 0.75rem;
        font-weight: 700;
        margin-bottom: 0.5rem;
    }
    .step-active { background: #22c55e; color: #fff; }
    .step-done { background: #3b82f6; color: #fff; }
    .step-locked { background: #374151; color: #9ca3af; }

    /* Sequence display */
    .seq-box {
        font-family: 'JetBrains Mono', monospace;
        font-size: 0.78rem;
        background: #0f172a;
        color: #a5f3fc;
        padding: 12px 16px;
        border-radius: 8px;
        word-break: break-all;
        line-height: 1.6;
        border: 1px solid #1e3a5f;
    }
    .domain-legend {
        display: inline-block;
        width: 14px; height: 14px;
        border-radius: 3px;
        margin-right: 6px;
        vertical-align: middle;
    }
</style>
""", unsafe_allow_html=True)

STEPS = [
    "1. Biomarker Input",
    "2. Nanobody Design",
    "3. Tag Selection",
    "4. Fusion Design",
    "5. CAR Build",
    "6. 3D Structure",
    "7. PK Simulation",
]


def step_status(idx: int) -> str:
    """Return CSS class for step status pill."""
    completed_keys = [
        "biomarker_info",
        "nanobody_candidate",
        "selected_tag",
        "fusion_construct",
        "car_construct",
        "structure_pdb",
        "pk_result",
    ]
    if idx < len(completed_keys) and completed_keys[idx] in st.session_state:
        return "step-done"
    current = st.session_state.get("current_step", 0)
    if idx == current:
        return "step-active"
    return "step-locked"


# ────────────────────────────────────────────
# Sidebar navigation
# ────────────────────────────────────────────
with st.sidebar:
    st.markdown("## 🧬 NanoCAR Designer")
    st.caption("Modular CAR-T Design Platform")
    st.markdown("---")

    if "current_step" not in st.session_state:
        st.session_state.current_step = 0

    for i, step_name in enumerate(STEPS):
        status = step_status(i)
        icon = "✅" if status == "step-done" else ("▶️" if status == "step-active" else "🔒")
        if st.button(f"{icon} {step_name}", key=f"nav_{i}", use_container_width=True):
            st.session_state.current_step = i

    st.markdown("---")
    st.markdown("##### Pipeline Progress")
    completed = sum(
        1 for k in [
            "biomarker_info", "nanobody_candidate", "selected_tag",
            "fusion_construct", "car_construct", "structure_pdb", "pk_result",
        ]
        if k in st.session_state
    )
    st.progress(completed / 7)
    st.caption(f"{completed}/7 steps completed")


# ────────────────────────────────────────────
# Helper functions
# ────────────────────────────────────────────
def show_sequence(seq: str, label: str = "Sequence") -> None:
    """Display a protein sequence in formatted box."""
    st.markdown(f"**{label}** ({len(seq)} aa)")
    st.markdown(f'<div class="seq-box">{seq}</div>', unsafe_allow_html=True)


def step_gate(required_key: str, required_step: str) -> bool:
    """Check if a previous step has been completed."""
    if required_key not in st.session_state:
        st.warning(f"⚠️ Please complete **{required_step}** first.")
        return False
    return True


# ════════════════════════════════════════════
# STEP 1: Biomarker Input
# ════════════════════════════════════════════
def render_step1() -> None:
    st.header("Step 1 — Biomarker Target Input")
    st.markdown(
        "Enter a **UniProt ID** to fetch protein data automatically, "
        "or paste a sequence manually."
    )

    tab_uniprot, tab_manual = st.tabs(["🔍 UniProt Lookup", "✏️ Manual Input"])

    with tab_uniprot:
        col1, col2 = st.columns([3, 1])
        with col1:
            uid = st.text_input(
                "UniProt Accession ID",
                placeholder="e.g. P00533 (EGFR), P04626 (HER2)",
                key="uniprot_input",
            )
        with col2:
            st.markdown("<br>", unsafe_allow_html=True)
            fetch_btn = st.button("Fetch", type="primary", key="fetch_uniprot")

        if fetch_btn and uid:
            with st.spinner("Fetching from UniProt..."):
                try:
                    info = fetch_from_uniprot(uid)
                    st.session_state.biomarker_info = info
                    st.success(f"✅ Retrieved: **{info.name}** ({info.gene_name})")
                except ValueError as e:
                    st.error(f"❌ {e}")
                except ConnectionError as e:
                    st.error(f"❌ Connection error: {e}")

    with tab_manual:
        name = st.text_input("Protein name", placeholder="e.g. CD19", key="manual_name")
        seq = st.text_area(
            "Amino acid sequence",
            height=120,
            placeholder="Paste single-letter amino acid sequence...",
            key="manual_seq",
        )
        if st.button("Validate & Use", key="manual_submit"):
            if seq:
                try:
                    info = parse_manual_input(seq, name or "User input")
                    st.session_state.biomarker_info = info
                    st.success(f"✅ Validated: {info.length} residues")
                except ValueError as e:
                    st.error(f"❌ {e}")

    # Display loaded biomarker info
    if "biomarker_info" in st.session_state:
        info = st.session_state.biomarker_info
        st.markdown("---")
        st.subheader(f"📋 {info.name}")

        c1, c2, c3 = st.columns(3)
        c1.metric("Length", f"{info.length} aa")
        c2.metric("Gene", info.gene_name or "N/A")
        c3.metric("Organism", info.organism or "N/A")

        if info.subcellular_location:
            st.info(f"📍 **Subcellular location:** {info.subcellular_location}")

        show_sequence(info.sequence, "Full Sequence")

        # Domains
        if info.domains:
            st.markdown("##### Domain Architecture")
            df = pd.DataFrame(info.domains)
            st.dataframe(df, use_container_width=True, hide_index=True)

        # Extracellular domains
        ec_domains = get_extracellular_domains(info)
        if ec_domains:
            st.markdown("##### 🎯 Extracellular Domains (targetable by nanobody)")
            for d in ec_domains:
                with st.expander(f"{d['description']} (pos {d['start']}-{d['end']}, {d['length']} aa)"):
                    show_sequence(d["sequence"], "Extracellular domain")

        # Disease associations
        if info.disease_associations:
            st.markdown("##### 🏥 Disease Associations")
            for disease in info.disease_associations:
                st.markdown(f"- {disease}")

        st.button(
            "Proceed to Step 2 →",
            on_click=lambda: setattr(st.session_state, "current_step", 1),
            type="primary",
        )


# ════════════════════════════════════════════
# STEP 2: Nanobody Design
# ════════════════════════════════════════════
def render_step2() -> None:
    st.header("Step 2 — Nanobody Design")

    if not step_gate("biomarker_info", "Step 1 (Biomarker Input)"):
        return

    info = st.session_state.biomarker_info
    st.markdown(f"**Target:** {info.name} ({info.gene_name})")

    tab_a, tab_b = st.tabs(["🔍 A: Database Search", "🧪 B: Template CDR Grafting"])

    # --- Tab A: DB search ---
    with tab_a:
        st.markdown("Search RCSB PDB for known VHH nanobodies against your target.")
        col1, col2 = st.columns([3, 1])
        with col1:
            search_term = st.text_input(
                "Search term",
                value=info.gene_name or info.name,
                key="nb_search_term",
            )
        with col2:
            st.markdown("<br>", unsafe_allow_html=True)
            search_btn = st.button("Search", type="primary", key="nb_search_btn")

        if search_btn and search_term:
            with st.spinner(f"Searching for VHH against '{search_term}'..."):
                candidates = search_known_nanobodies(search_term)
                st.session_state.nb_search_results = candidates

        if "nb_search_results" in st.session_state:
            results = st.session_state.nb_search_results
            if results:
                st.success(f"Found {len(results)} candidate(s)")
                for i, cand in enumerate(results):
                    with st.expander(f"Candidate {i+1} — PDB: {cand.pdb_id} ({len(cand.sequence)} aa)"):
                        show_sequence(cand.sequence, "VHH Sequence")
                        val = validate_nanobody(cand.sequence)
                        for p in val["passed"]:
                            st.markdown(f"✅ {p}")
                        for w in val["warnings"]:
                            st.markdown(f"⚠️ {w}")
                        for e in val["errors"]:
                            st.markdown(f"❌ {e}")

                        if st.button(f"Select Candidate {i+1}", key=f"select_nb_{i}"):
                            cand.validation = val
                            cand.target_antigen = search_term
                            regions = annotate_regions(cand.sequence)
                            if "error" not in regions:
                                cand.cdr1 = regions["CDR1"]["sequence"]
                                cand.cdr2 = regions["CDR2"]["sequence"]
                                cand.cdr3 = regions["CDR3"]["sequence"]
                            st.session_state.nanobody_candidate = cand
                            st.rerun()
            else:
                st.info("No VHH candidates found. Try Tab B (Template CDR Grafting).")

    # --- Tab B: Template CDR grafting ---
    with tab_b:
        st.markdown("Graft custom CDR sequences onto a humanized VHH framework.")

        templates = load_framework_templates()
        guidelines = get_cdr_guidelines()

        fw_options = {t["id"]: f"{t['name']} — {t['description'][:60]}..." for t in templates}
        selected_fw = st.selectbox(
            "Framework template",
            options=list(fw_options.keys()),
            format_func=lambda x: fw_options[x],
            key="fw_select",
        )

        # Show framework details
        fw = next(t for t in templates if t["id"] == selected_fw)
        with st.expander("Framework details"):
            st.markdown(f"**Organism:** {fw['organism_origin']}")
            st.markdown(f"**Tm:** {fw['thermal_stability_tm']}°C")
            st.markdown(f"**Expression hosts:** {', '.join(fw['expression_host'])}")
            st.code(f"FR1: {fw['fr1']}\nFR2: {fw['fr2']}\nFR3: {fw['fr3']}\nFR4: {fw['fr4']}")

        st.markdown("##### CDR Sequences")
        st.caption(
            f"CDR1: {guidelines['CDR1']['typical_length']} aa | "
            f"CDR2: {guidelines['CDR2']['typical_length']} aa | "
            f"CDR3: {guidelines['CDR3']['typical_length']} aa (can be longer in VHH)"
        )

        cdr1 = st.text_input("CDR1", placeholder="e.g. GFTFSSYA", key="cdr1_input")
        cdr2 = st.text_input("CDR2", placeholder="e.g. ISGSGGST", key="cdr2_input")
        cdr3 = st.text_input("CDR3", placeholder="e.g. AKDRGYRYYGSSHWYFDV", key="cdr3_input")

        if st.button("Build Nanobody", type="primary", key="build_nb"):
            if cdr1 and cdr2 and cdr3:
                try:
                    candidate = graft_cdrs(selected_fw, cdr1, cdr2, cdr3)
                    candidate.target_antigen = info.gene_name or info.name
                    st.session_state.nanobody_candidate = candidate
                    st.success("✅ Nanobody assembled successfully")
                    st.rerun()
                except ValueError as e:
                    st.error(f"❌ {e}")
            else:
                st.error("Please enter all three CDR sequences")

    # Display selected nanobody
    if "nanobody_candidate" in st.session_state:
        nb = st.session_state.nanobody_candidate
        st.markdown("---")
        st.subheader("✅ Selected Nanobody")

        c1, c2, c3 = st.columns(3)
        c1.metric("Source", nb.source.title())
        c2.metric("Length", f"{len(nb.sequence)} aa")
        c3.metric("Target", nb.target_antigen)

        show_sequence(nb.sequence, "Full VHH Sequence")

        if nb.cdr1:
            st.markdown(
                f"**CDR1:** `{nb.cdr1}` | **CDR2:** `{nb.cdr2}` | **CDR3:** `{nb.cdr3}`"
            )

        st.button(
            "Proceed to Step 3 →",
            on_click=lambda: setattr(st.session_state, "current_step", 2),
            type="primary",
        )


# ════════════════════════════════════════════
# STEP 3: Tag Selection
# ════════════════════════════════════════════
def render_step3() -> None:
    st.header("Step 3 — Tag System Selection")

    tags = load_tag_database()
    tag_names = [t.name for t in tags]

    # Radar chart comparison
    st.markdown("##### Comparative Radar Chart")
    selected_for_chart = st.multiselect(
        "Select tags to compare",
        tag_names,
        default=tag_names,
        key="radar_tags",
    )

    if selected_for_chart:
        radar_data = get_radar_chart_data(selected_for_chart)
        fig = go.Figure()
        colors = ["#22c55e", "#3b82f6", "#f59e0b", "#ef4444", "#8b5cf6"]
        for i, trace in enumerate(radar_data["traces"]):
            fig.add_trace(go.Scatterpolar(
                r=trace["values"],
                theta=trace["categories"],
                fill="toself",
                name=trace["name"],
                line_color=colors[i % len(colors)],
                opacity=0.7,
            ))
        fig.update_layout(
            polar=dict(
                radialaxis=dict(visible=True, range=[0, 10]),
                bgcolor="rgba(0,0,0,0)",
            ),
            showlegend=True,
            height=420,
            margin=dict(l=60, r=60, t=30, b=30),
            paper_bgcolor="rgba(0,0,0,0)",
            plot_bgcolor="rgba(0,0,0,0)",
        )
        st.plotly_chart(fig, use_container_width=True)

    # Recommendation
    st.markdown("##### 🎯 Quick Recommendation")
    priority = st.selectbox(
        "Priority",
        ["safety", "efficacy", "clinical_readiness", "manufacturability"],
        key="tag_priority",
    )
    rec = recommend_tag(priority)
    st.info(f"Recommended tag for **{priority}**: **{rec.name}** ({rec.platform})")

    # Comparison table
    st.markdown("##### Detailed Comparison")
    comp = compare_tags(tag_names)
    df_comp = pd.DataFrame(comp["comparison_table"])
    st.dataframe(df_comp, use_container_width=True, hide_index=True)

    # Individual tag details
    st.markdown("##### Select a Tag System")
    selected_tag_name = st.selectbox(
        "Choose tag",
        tag_names,
        format_func=lambda x: f"{x} — {get_tag_by_name(x).platform}",
        key="tag_select",
    )

    tag = get_tag_by_name(selected_tag_name)
    if tag:
        col1, col2 = st.columns(2)
        with col1:
            st.markdown(f"**Full name:** {tag.full_name}")
            st.markdown(f"**Platform:** {tag.platform}")
            st.markdown(f"**Origin:** {tag.origin}")
            st.markdown(f"**MW:** {tag.molecular_weight_da} Da")
            st.markdown(f"**Half-life:** {tag.pk_half_life_min} min")
            st.markdown(f"**Immunogenicity:** {tag.immunogenicity}")
            st.markdown(f"**Clinical stage:** {tag.clinical_stage}")
            if tag.sequence != "N/A (small molecule conjugate)":
                show_sequence(tag.sequence, "Tag Sequence")
        with col2:
            st.markdown("**✅ Pros:**")
            for p in tag.pros:
                st.markdown(f"- {p}")
            st.markdown("**⚠️ Cons:**")
            for c in tag.cons:
                st.markdown(f"- {c}")

        if st.button("Select This Tag", type="primary", key="confirm_tag"):
            st.session_state.selected_tag = tag
            st.success(f"✅ Selected: {tag.name} ({tag.platform})")
            st.rerun()

    if "selected_tag" in st.session_state:
        st.markdown("---")
        sel = st.session_state.selected_tag
        st.success(f"**Current selection:** {sel.name} ({sel.platform})")
        st.button(
            "Proceed to Step 4 →",
            on_click=lambda: setattr(st.session_state, "current_step", 3),
            type="primary",
        )


# ════════════════════════════════════════════
# STEP 4: Fusion Design
# ════════════════════════════════════════════
def render_step4() -> None:
    st.header("Step 4 — Nb-Tag Fusion Construct Design")

    if not step_gate("nanobody_candidate", "Step 2 (Nanobody Design)"):
        return
    if not step_gate("selected_tag", "Step 3 (Tag Selection)"):
        return

    nb = st.session_state.nanobody_candidate
    tag = st.session_state.selected_tag

    st.markdown(f"**Nanobody:** {len(nb.sequence)} aa | **Tag:** {tag.name} ({tag.platform})")

    # Design parameters
    col1, col2 = st.columns(2)
    with col1:
        linkers = load_linker_library()
        linker_options = {ln["id"]: f"{ln['name']} — {ln['type']} ({ln['length_aa']} aa)" for ln in linkers}
        selected_linker = st.selectbox(
            "Linker",
            options=list(linker_options.keys()),
            format_func=lambda x: linker_options[x],
            index=2,  # default to (G4S)x3
            key="linker_select",
        )

        # Show linker details
        ln = next(l for l in linkers if l["id"] == selected_linker)
        st.caption(f"Type: {ln['type']} | ~{ln['approx_length_nm']} nm | Seq: `{ln['sequence']}`")

    with col2:
        orientation = st.selectbox(
            "Orientation",
            ["Nb-linker-Tag", "Tag-linker-Nb"],
            key="orient_select",
        )
        codon_org = st.selectbox(
            "Codon optimization",
            ["human", "ecoli"],
            key="codon_select",
        )

    if st.button("Design Fusion Construct", type="primary", key="design_fusion"):
        tag_seq = tag.sequence if tag.sequence != "N/A (small molecule conjugate)" else "FITC"
        try:
            fusion = design_fusion(
                nanobody_seq=nb.sequence,
                tag_seq=tag_seq,
                linker_id=selected_linker,
                orientation=orientation,
                codon_organism=codon_org,
            )
            st.session_state.fusion_construct = fusion
            st.success("✅ Fusion construct designed")
            st.rerun()
        except ValueError as e:
            st.error(f"❌ {e}")

    if "fusion_construct" in st.session_state:
        fusion = st.session_state.fusion_construct
        st.markdown("---")
        st.subheader("📐 Fusion Construct")

        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Total length", f"{fusion.total_length} aa")
        c2.metric("MW", f"{fusion.molecular_weight/1000:.1f} kDa")
        c3.metric("pI", f"{fusion.isoelectric_point:.2f}")
        c4.metric("GRAVY", f"{fusion.gravy:.3f}")

        stability = "Stable" if fusion.instability_index < 40 else "Potentially unstable"
        st.markdown(f"**Instability Index:** {fusion.instability_index:.1f} ({stability})")

        show_sequence(fusion.full_sequence, f"Fusion Protein ({fusion.orientation})")

        # Domain map visualization
        nb_len = len(fusion.nanobody_seq)
        lk_len = len(fusion.linker_seq)
        tag_len = len(fusion.tag_seq)

        st.markdown("##### Domain Map")
        domain_bar = ""
        if fusion.orientation == "Nb-linker-Tag":
            domain_bar = (
                f'<span class="domain-legend" style="background:#2ecc71"></span>Nanobody ({nb_len} aa) | '
                f'<span class="domain-legend" style="background:#f39c12"></span>Linker ({lk_len} aa) | '
                f'<span class="domain-legend" style="background:#e74c3c"></span>Tag ({tag_len} aa)'
            )
        else:
            domain_bar = (
                f'<span class="domain-legend" style="background:#e74c3c"></span>Tag ({tag_len} aa) | '
                f'<span class="domain-legend" style="background:#f39c12"></span>Linker ({lk_len} aa) | '
                f'<span class="domain-legend" style="background:#2ecc71"></span>Nanobody ({nb_len} aa)'
            )
        st.markdown(domain_bar, unsafe_allow_html=True)

        # DNA sequence
        with st.expander("🧬 DNA sequence (codon-optimized)"):
            st.code(fusion.dna_sequence, language=None)
            st.caption(f"Length: {len(fusion.dna_sequence)} bp")

        st.button(
            "Proceed to Step 5 →",
            on_click=lambda: setattr(st.session_state, "current_step", 4),
            type="primary",
        )


# ════════════════════════════════════════════
# STEP 5: CAR Build
# ════════════════════════════════════════════
def render_step5() -> None:
    st.header("Step 5 — CAR Construct Design")

    if not step_gate("selected_tag", "Step 3 (Tag Selection)"):
        return

    tag = st.session_state.selected_tag
    st.markdown(f"**Tag system:** {tag.name} ({tag.platform}) → Recognition: {tag.recognition_domain_name}")

    available = list_available_domains()
    lib = load_domain_library()

    col1, col2 = st.columns(2)
    with col1:
        hinge_key = st.selectbox(
            "Hinge domain",
            available["hinge_domains"],
            format_func=lambda x: lib["hinge_domains"][x]["name"],
            key="hinge_sel",
        )
        tm_key = st.selectbox(
            "Transmembrane domain",
            available["transmembrane_domains"],
            format_func=lambda x: lib["transmembrane_domains"][x]["name"],
            key="tm_sel",
        )
    with col2:
        costim_key = st.selectbox(
            "Costimulatory domain",
            available["costimulatory_domains"],
            format_func=lambda x: lib["costimulatory_domains"][x]["name"],
            key="costim_sel",
        )
        sig_key = st.selectbox(
            "Signaling domain",
            available["signaling_domains"],
            format_func=lambda x: lib["signaling_domains"][x]["name"],
            key="sig_sel",
        )

    if st.button("Build CAR Construct", type="primary", key="build_car"):
        try:
            car = build_car(
                tag_name=tag.name,
                hinge_key=hinge_key,
                tm_key=tm_key,
                costim_key=costim_key,
                signaling_key=sig_key,
            )
            st.session_state.car_construct = car
            st.success("✅ CAR construct assembled")
            st.rerun()
        except ValueError as e:
            st.error(f"❌ {e}")

    if "car_construct" in st.session_state:
        car = st.session_state.car_construct
        summary = get_car_summary(car)
        st.markdown("---")
        st.subheader("🔬 CAR Construct")

        c1, c2, c3 = st.columns(3)
        c1.metric("Total length", f"{summary['total_length_aa']} aa")
        c2.metric("Generation", summary["generation"])
        c3.metric("Costimulatory", summary["costimulatory"])

        # Domain table
        st.markdown("##### Domain Architecture")
        df = pd.DataFrame(summary["domains"])
        st.dataframe(df, use_container_width=True, hide_index=True)

        # Domain map visual
        st.markdown("##### Domain Map")
        colors = ["#2ecc71", "#3498db", "#9b59b6", "#f39c12", "#e74c3c"]
        legend = ""
        for i, d in enumerate(car.domain_map):
            legend += (
                f'<span class="domain-legend" style="background:{colors[i]}"></span>'
                f'{d["domain"]} ({d["length"]} aa) | '
            )
        st.markdown(legend.rstrip(" | "), unsafe_allow_html=True)

        show_sequence(car.full_sequence, "Full CAR Sequence")

        st.button(
            "Proceed to Step 6 →",
            on_click=lambda: setattr(st.session_state, "current_step", 5),
            type="primary",
        )


# ════════════════════════════════════════════
# STEP 6: 3D Structure
# ════════════════════════════════════════════
def render_step6() -> None:
    st.header("Step 6 — Structure Prediction & 3D Visualization")

    has_fusion = "fusion_construct" in st.session_state
    has_car = "car_construct" in st.session_state

    if not has_fusion and not has_car:
        st.warning("⚠️ Please complete Step 4 or Step 5 first.")
        return

    target = st.selectbox(
        "Predict structure for",
        [x for x in ["Fusion Construct", "CAR Construct"] if
         (x == "Fusion Construct" and has_fusion) or
         (x == "CAR Construct" and has_car)],
        key="struct_target",
    )

    if target == "Fusion Construct":
        seq = st.session_state.fusion_construct.full_sequence
    else:
        seq = st.session_state.car_construct.full_sequence

    st.caption(f"Sequence length: {len(seq)} aa")
    if len(seq) > 800:
        st.warning("⚠️ Sequence exceeds 800 aa — will be truncated for prediction.")

    if st.button("🚀 Predict Structure", type="primary", key="predict_struct"):
        with st.spinner("Running structure prediction (ESMFold)... This may take 1-2 minutes."):
            pdb = predict_structure(seq)
            if pdb:
                st.session_state.structure_pdb = {
                    "pdb_data": pdb,
                    "target": target,
                }
                st.success("✅ Structure predicted successfully")
                st.rerun()
            else:
                st.error("❌ Structure prediction failed. The API may be temporarily unavailable.")

    if "structure_pdb" in st.session_state:
        pdb_data = st.session_state.structure_pdb["pdb_data"]
        target_type = st.session_state.structure_pdb["target"]

        st.markdown("---")
        st.subheader(f"🧊 3D Structure — {target_type}")

        color_scheme = st.selectbox(
            "Color scheme",
            ["spectrum", "domain", "hydrophobicity"],
            key="color_scheme",
        )

        if color_scheme == "domain" and target_type == "Fusion Construct":
            f = st.session_state.fusion_construct
            html = visualize_fusion_construct(
                pdb_data,
                nb_length=len(f.nanobody_seq),
                linker_length=len(f.linker_seq),
                tag_length=len(f.tag_seq),
            )
        elif color_scheme == "domain" and target_type == "CAR Construct":
            car = st.session_state.car_construct
            html = visualize_car_domains(pdb_data, car.domain_map)
        else:
            html = render_3d_viewer(pdb_data, color_scheme=color_scheme)

        components.html(html, height=520)

        # Download PDB
        st.download_button(
            "📥 Download PDB",
            data=pdb_data,
            file_name=f"nanocar_{target_type.lower().replace(' ', '_')}.pdb",
            mime="text/plain",
        )

        st.button(
            "Proceed to Step 7 →",
            on_click=lambda: setattr(st.session_state, "current_step", 6),
            type="primary",
        )


# ════════════════════════════════════════════
# STEP 7: PK Simulation
# ════════════════════════════════════════════
def render_step7() -> None:
    st.header("Step 7 — PK / Safety Simulation")

    tag = st.session_state.get("selected_tag")
    default_hl = tag.pk_half_life_min if tag else 15.0

    st.markdown("##### Pharmacokinetic Parameters")
    col1, col2 = st.columns(2)
    with col1:
        half_life = st.number_input("Half-life (min)", 1.0, 600.0, default_hl, key="pk_hl")
        dose = st.number_input("Dose (mg/kg)", 0.01, 100.0, 1.0, key="pk_dose")
        interval = st.number_input("Dosing interval (hr)", 0.5, 48.0, 4.0, key="pk_interval")
        num_doses = st.number_input("Number of doses", 1, 20, 6, key="pk_ndose")
    with col2:
        vd = st.number_input("Volume of distribution (L)", 0.5, 50.0, 5.0, key="pk_vd")
        body_wt = st.number_input("Body weight (kg)", 20.0, 150.0, 70.0, key="pk_bw")
        ec50 = st.number_input("EC50 (ng/mL)", 1.0, 5000.0, 100.0, key="pk_ec50")
        safety_th = st.number_input("Safety threshold (ng/mL)", 10.0, 50000.0, 1000.0, key="pk_safety")

    if st.button("🚀 Run Simulation", type="primary", key="run_pk"):
        params = PKParameters(
            half_life_min=half_life,
            dose_mg_kg=dose,
            dosing_interval_hr=interval,
            num_doses=int(num_doses),
            volume_distribution_L=vd,
            body_weight_kg=body_wt,
        )

        with st.spinner("Running PK simulation..."):
            result = simulate_on_off_dynamics(params, ec50, safety_th)
            st.session_state.pk_result = result
            st.rerun()

    if "pk_result" in st.session_state:
        result = st.session_state.pk_result
        pk = result["pk_result"]
        report = generate_pk_report(result)

        st.markdown("---")
        st.subheader("📊 Simulation Results")

        # Key metrics
        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Cmax", f"{report['Cmax (ng/mL)']} ng/mL")
        c2.metric("Trough", f"{report['Trough (ng/mL)']} ng/mL")
        c3.metric("ON time", report["ON time fraction"])
        c4.metric("Time to OFF", f"{report['Time to OFF after last dose (min)']} min")

        # Concentration-time plot
        st.markdown("##### Concentration-Time Profile")
        time_hrs = pk.time_points
        conc = pk.concentration

        # Downsample for plotly (every 10th point)
        step = max(1, len(time_hrs) // 500)
        t_ds = time_hrs[::step]
        c_ds = conc[::step]
        a_ds = result["activation"][::step]

        fig = go.Figure()
        fig.add_trace(go.Scatter(
            x=t_ds, y=c_ds,
            mode="lines",
            name="Concentration",
            line=dict(color="#3b82f6", width=2),
        ))

        # EC50 and safety lines
        fig.add_hline(y=result["ec50"], line_dash="dash",
                      line_color="#22c55e", annotation_text="EC50")
        fig.add_hline(y=result["safety_threshold"], line_dash="dash",
                      line_color="#ef4444", annotation_text="Safety limit")

        fig.update_layout(
            xaxis_title="Time (hours)",
            yaxis_title="Concentration (ng/mL)",
            height=380,
            margin=dict(l=50, r=30, t=30, b=50),
            paper_bgcolor="rgba(0,0,0,0)",
            plot_bgcolor="rgba(0,0,0,0)",
        )
        st.plotly_chart(fig, use_container_width=True)

        # Activation dynamics
        st.markdown("##### ON/OFF Activation Dynamics")
        fig2 = go.Figure()
        fig2.add_trace(go.Scatter(
            x=t_ds, y=[a * 100 for a in a_ds],
            mode="lines",
            name="CAR-T Activation",
            fill="tozeroy",
            line=dict(color="#22c55e", width=2),
            fillcolor="rgba(34,197,94,0.2)",
        ))
        fig2.add_hline(y=50, line_dash="dot",
                       line_color="#f59e0b", annotation_text="50% activation")
        fig2.update_layout(
            xaxis_title="Time (hours)",
            yaxis_title="Activation (%)",
            yaxis=dict(range=[0, 105]),
            height=300,
            margin=dict(l=50, r=30, t=30, b=50),
            paper_bgcolor="rgba(0,0,0,0)",
            plot_bgcolor="rgba(0,0,0,0)",
        )
        st.plotly_chart(fig2, use_container_width=True)

        # Safety assessment
        st.markdown("##### Safety Assessment")
        violation = result["safety_violation_fraction"]
        if violation == 0:
            st.success("✅ No safety threshold violations detected.")
        elif violation < 0.05:
            st.warning(f"⚠️ Minor safety violations ({violation*100:.1f}% of time above threshold)")
        else:
            st.error(f"❌ Significant safety violations ({violation*100:.1f}% of time above threshold)")

        # Full report table
        with st.expander("📋 Full PK Report"):
            df_report = pd.DataFrame(
                [{"Parameter": k, "Value": v} for k, v in report.items()]
            )
            st.dataframe(df_report, use_container_width=True, hide_index=True)

        st.success("🎉 **Pipeline Complete!** All 7 steps have been executed.")


# ════════════════════════════════════════════
# Main routing
# ════════════════════════════════════════════
STEP_RENDERERS = [
    render_step1,
    render_step2,
    render_step3,
    render_step4,
    render_step5,
    render_step6,
    render_step7,
]

current = st.session_state.get("current_step", 0)
STEP_RENDERERS[current]()
