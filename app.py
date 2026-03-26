#!/usr/bin/env python3
"""
SAIDock Web Interface
Streamlit frontend that communicates with FastAPI backend.
For local testing without the backend, uses direct pipeline call.
"""
import streamlit as st
import requests
import json
import time
import os
from pathlib import Path

API_URL = os.environ.get('SAIDOCK_API_URL', 'http://localhost:8000')

st.set_page_config(
    page_title = 'SAIDock',
    page_icon  = 'S',
    layout     = 'wide',
    initial_sidebar_state = 'expanded',
)

st.markdown("""
<style>
  body       { font-family: Arial, sans-serif; }
  .header    { background: #1B2A4A; color: white;
               padding: 20px 30px; border-radius: 8px;
               margin-bottom: 20px; }
  .header h1 { margin: 0; font-size: 1.8em; }
  .header p  { margin: 5px 0 0 0; opacity: 0.85; }
  .dtss-card { border-radius: 8px; padding: 14px;
               text-align: center; color: white;
               font-size: 2.0em; font-weight: bold; }
  .metric-card { background: #F8F9FA; border-radius: 6px;
                 padding: 12px; text-align: center;
                 border: 1px solid #E0E0E0; }
</style>
""", unsafe_allow_html=True)

st.markdown("""
<div class="header">
  <h1>SAIDock v1.0.0</h1>
  <p>Automated Drug-Target Docking and Assessment Pipeline |
     Pocket Detection | ADMET | ML Scoring | DTSS</p>
</div>
""", unsafe_allow_html=True)

# ---- Sidebar ----------------------------------------------------------------
with st.sidebar:
    st.markdown("### Ligand Input")
    input_mode = st.radio("Input type",
                          ["Drug name", "SMILES", "Upload SDF"],
                          index=0)
    ligand_val = None
    if input_mode == "Drug name":
        ligand_val = st.text_input("Compound name",
                                   placeholder="Ellagic acid, Quercetin, ...")
        st.caption("Fetches SMILES automatically from PubChem")
    elif input_mode == "SMILES":
        ligand_val = st.text_area("SMILES", height=90,
                                  placeholder="COc1cc(/C=C/C(=O)O)ccc1O")
    else:
        sdf_file = st.file_uploader("Upload SDF", type=['sdf','mol'])
        if sdf_file:
            p = Path(f'/tmp/{sdf_file.name}')
            p.write_bytes(sdf_file.read())
            ligand_val = str(p)

    st.markdown("### Target Input")
    target_mode = st.radio("Target type", ["PDB ID", "Upload PDB"])
    if target_mode == "PDB ID":
        target_val = st.text_input("PDB ID",
                                   placeholder="3PE1, 5L2I, 4XMB ...")
        st.caption("Downloads from RCSB automatically")
    else:
        pdb_file = st.file_uploader("Upload PDB", type=['pdb'])
        target_val = None
        if pdb_file:
            p = Path(f'/tmp/{pdb_file.name}')
            p.write_bytes(pdb_file.read())
            target_val = str(p)

    st.markdown("### Parameters")
    with st.expander("Docking settings"):
        exhaustiveness = st.slider("Exhaustiveness", 4, 64, 16)
        n_pockets      = st.slider("Pockets to screen", 1, 8, 3)
        chain          = st.text_input("Receptor chain", "A")
        cpu            = st.slider("CPU cores", 1, 16, 4)
    with st.expander("Analysis options"):
        do_admet = st.checkbox("ADMET assessment",       value=True)
        do_inter = st.checkbox("Interaction analysis",   value=True)
        do_ml    = st.checkbox("ML scoring and DTSS",    value=True)

    run_btn = st.button("Run SAIDock Assessment",
                        type="primary", use_container_width=True)

# ---- Main area -------------------------------------------------------------
if not run_btn:
    col1, col2, col3, col4 = st.columns(4)
    panels = [
        ("Target Preparation",
         "PDB download, chain extraction, missing residue repair, PDBQT"),
        ("Pocket Detection",
         "fpocket or PolarPocket: alpha sphere detection, druggability ranking"),
        ("Docking",
         "AutoDock Vina or PolarDock: multi-pocket screening, pose clustering"),
        ("DTSS Score",
         "Novel composite: DG + ADMET + ML + Pocket druggability + Literature"),
    ]
    for col, (title, desc) in zip([col1,col2,col3,col4], panels):
        with col:
            st.markdown(f"""
            <div class="metric-card">
              <b>{title}</b><br>
              <small style="color:#555">{desc}</small>
            </div>
            """, unsafe_allow_html=True)

    st.markdown("---")
    st.markdown("### Quick Examples")
    c1, c2, c3 = st.columns(3)
    with c1:
        if st.button("Ferulic Acid vs CK2a (3PE1)",
                     use_container_width=True):
            st.session_state['ql'] = 'COc1cc(/C=C/C(=O)O)ccc1O'
            st.session_state['qt'] = '3PE1'
    with c2:
        if st.button("Caffeic Acid vs Keap1 (4XMB)",
                     use_container_width=True):
            st.session_state['ql'] = 'OC(=O)/C=C/c1ccc(O)c(O)c1'
            st.session_state['qt'] = '4XMB'
    with c3:
        if st.button("Ellagic Acid vs CDK6 (5L2I)",
                     use_container_width=True):
            st.session_state['ql'] = \
              'OC1=C2C(=CC(=C1)O)C1=C(OC2=O)C(=O)OC2=CC(=CC(=C21)O)O'
            st.session_state['qt'] = '5L2I'

else:
    if not ligand_val:
        st.error("Please provide a ligand input in the sidebar.")
        st.stop()
    if not target_val:
        st.error("Please provide a target (PDB ID or file) in the sidebar.")
        st.stop()

    st.markdown("---")
    progress_bar = st.progress(0)
    status_msg   = st.empty()
    log_area     = st.empty()

    payload = {
        'ligand_name':    ligand_val if input_mode == "Drug name" else None,
        'smiles':         ligand_val if input_mode == "SMILES"    else None,
        'sdf_path':       ligand_val if input_mode == "Upload SDF" else None,
        'target':         target_val,
        'chain':          chain,
        'exhaustiveness': exhaustiveness,
        'n_pockets':      n_pockets,
        'cpu':            cpu,
        'do_admet':       do_admet,
        'do_interactions':do_inter,
        'do_ml':          do_ml,
    }

    # Try API first, fall back to direct local run
    try:
        status_msg.info("Submitting job to SAIDock server...")
        resp   = requests.post(f'{API_URL}/submit',
                               json=payload, timeout=10)
        job    = resp.json()
        job_id = job['job_id']
        status_msg.success(f"Job submitted. ID: {job_id}")
        _run_mode = 'api'
    except Exception:
        status_msg.warning(
            "API server not reachable. Running directly in this session."
        )
        _run_mode = 'local'

    if _run_mode == 'api':
        while True:
            try:
                stat = requests.get(f'{API_URL}/status/{job_id}',
                                    timeout=10).json()
            except Exception:
                time.sleep(5)
                continue
            pct  = stat.get('progress_pct', 0)
            step = stat.get('current_step', '...')
            progress_bar.progress(int(pct))
            status_msg.info(f"[{pct}%] {step}")
            lines = stat.get('log_lines', [])
            if lines:
                log_area.code('\n'.join(lines[-15:]))
            if stat['status'] == 'COMPLETE':
                progress_bar.progress(100)
                status_msg.success("Assessment complete.")
                _show_results(stat['result'])
                st.link_button(
                    "Download Full HTML Report",
                    f"{API_URL}/download/{job_id}/report",
                    use_container_width=True,
                )
                break
            elif stat['status'] == 'FAILED':
                status_msg.error(f"Failed: {stat.get('error','')}")
                break
            time.sleep(4)

    else:
        import types, sys
        sys.path.insert(0, str(Path(__file__).parent.parent))
        from saidock.pipeline import SAIDockRun
        args = types.SimpleNamespace(
            ligand          = payload.get('ligand_name'),
            smiles          = payload.get('smiles'),
            sdf             = payload.get('sdf_path'),
            target          = payload['target'],
            chain           = payload.get('chain','A'),
            output          = '/tmp/saidock_streamlit_run/',
            exhaustiveness  = payload.get('exhaustiveness',16),
            n_pockets       = payload.get('n_pockets',3),
            cpu             = payload.get('cpu',4),
            no_admet        = not payload.get('do_admet',True),
            no_interactions = not payload.get('do_interactions',True),
            no_ml           = not payload.get('do_ml',True),
        )
        progress_bar.progress(10)
        status_msg.info("Running SAIDock pipeline...")
        try:
            runner = SAIDockRun(args)
            runner.execute()
            progress_bar.progress(100)
            status_msg.success("Assessment complete.")
            _show_results(runner.state)
            report_path = Path('/tmp/saidock_streamlit_run/report/saidock_report.html')
            if report_path.exists():
                st.download_button(
                    label    = "Download HTML Report",
                    data     = report_path.read_bytes(),
                    file_name= 'saidock_report.html',
                    mime     = 'text/html',
                    use_container_width=True,
                )
        except Exception as e:
            status_msg.error(f"Pipeline error: {e}")
            import traceback
            st.code(traceback.format_exc())


def _show_results(state: dict):
    ml  = state.get('ml_results',  {})
    amd = state.get('admet',       {})
    dr  = state.get('docking_results', [{}])
    dtss = ml.get('DTSS', 0)
    cat  = ml.get('binding_category', 'N/A')
    best = dr[0].get('best_score', 0) if dr else 0

    color_map = {
        'Excellent':'#27AE60','Good':'#2980B9',
        'Moderate':'#F39C12','Weak':'#E74C3C','Poor':'#95A5A6',
    }
    color = color_map.get(cat.split()[0] if cat else 'Poor', '#95A5A6')

    st.markdown("---")
    st.markdown("## Assessment Results")

    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.markdown(f"""
        <div class="dtss-card" style="background:{color}">
          DTSS = {dtss:.3f}
        </div>
        <div style="text-align:center; margin-top:6px">
          <small>{cat}</small>
        </div>
        """, unsafe_allow_html=True)
    with col2:
        st.metric("Best Vina Score",
                  f"{best:.2f} kcal/mol")
    with col3:
        st.metric("QED", f"{amd.get('QED', 0):.3f}")
    with col4:
        st.metric("Pockets Screened", len(dr))

    import plotly.graph_objects as go
    import pandas as pd

    comps = {
        'DG (x0.35)':         ml.get('dg_norm',0)*0.35,
        'ADMET (x0.20)':      ml.get('admet_score',0)*0.20,
        'ML conf (x0.20)':    ml.get('ml_confidence',0)*0.20,
        'Pocket (x0.15)':     ml.get('pocket_druggability',0)*0.15,
        'Literature (x0.10)': ml.get('literature_score',0)*0.10,
    }
    fig = go.Figure(go.Bar(
        x=list(comps.values()), y=list(comps.keys()),
        orientation='h',
        text=[f'{v:.3f}' for v in comps.values()],
        textposition='outside',
        marker_color=['#1B4F8A','#1A7A4A','#9B59B6','#F39C12','#E74C3C'],
    ))
    fig.add_vline(x=dtss, line_dash='dash', line_color='#C0392B',
                  annotation_text=f'DTSS={dtss:.3f}')
    fig.update_layout(
        title='DTSS Component Breakdown',
        height=240,
        margin=dict(l=150,r=90,t=40,b=30),
        plot_bgcolor='white', paper_bgcolor='white',
    )
    st.plotly_chart(fig, use_container_width=True)

    col_d, col_a = st.columns(2)
    with col_d:
        st.markdown("#### Docking Results by Pocket")
        dock_df = pd.DataFrame([{
            'Pocket': f'Pocket {r.get("pocket_id","?")}',
            'Score (kcal/mol)': r.get('best_score', 0),
            'Mean Top-3':       r.get('mean_top3', 0),
            'Engine':           r.get('engine', ''),
        } for r in dr])
        st.dataframe(dock_df, use_container_width=True)
    with col_a:
        st.markdown("#### ADMET Properties")
        admet_df = pd.DataFrame([
            {'Property': 'MW',       'Value': amd.get('MW',''),
             'Limit': 'le 500',
             'Status': 'Pass' if float(amd.get('MW',0) or 0) <= 500 else 'Fail'},
            {'Property': 'LogP',     'Value': amd.get('LogP',''),
             'Limit': 'le 5',
             'Status': 'Pass' if float(amd.get('LogP',0) or 0) <= 5 else 'Fail'},
            {'Property': 'TPSA',     'Value': amd.get('TPSA',''),
             'Limit': 'le 140',
             'Status': 'Pass' if float(amd.get('TPSA',0) or 0) <= 140 else 'Fail'},
            {'Property': 'QED',      'Value': amd.get('QED',''),
             'Limit': 'gt 0.5', 'Status': ''},
            {'Property': 'GI Abs',   'Value': amd.get('GI_absorption',''),
             'Limit': 'High',   'Status': ''},
            {'Property': 'BBB',      'Value': str(amd.get('BBB_penetrant','')),
             'Limit': '',        'Status': ''},
            {'Property': 'Solubility','Value': amd.get('Solubility_class',''),
             'Limit': '',         'Status': ''},
        ])
        st.dataframe(admet_df, use_container_width=True)
