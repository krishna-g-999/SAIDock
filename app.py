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


# ── SAIDock Header (matches local server exactly) ────────────────────
st.markdown("""
<style>
    /* Main header bar */
    .saidock-header {
        background: linear-gradient(135deg, #0D2137 0%, #1a3a5c 100%);
        padding: 2rem 2.5rem 1.5rem 2.5rem;
        border-radius: 8px;
        margin-bottom: 1.5rem;
    }
    .saidock-header h1 {
        color: #FFFFFF;
        font-size: 2rem;
        font-weight: 700;
        letter-spacing: 0.02em;
        margin: 0 0 0.4rem 0;
        font-family: 'Source Sans Pro', sans-serif;
    }
    .saidock-header p {
        color: #A8C4E0;
        font-size: 0.95rem;
        margin: 0;
        font-family: 'Source Sans Pro', sans-serif;
    }
    /* Sidebar styling */
    [data-testid="stSidebar"] {
        background-color: #F0F4F8;
    }
    [data-testid="stSidebar"] .stRadio label {
        font-size: 0.9rem;
        color: #0D2137;
    }
    /* Run button */
    .stButton > button {
        background-color: #0D2137;
        color: white;
        border: none;
        border-radius: 6px;
        font-weight: 600;
        width: 100%;
        padding: 0.6rem 1rem;
    }
    .stButton > button:hover {
        background-color: #1a3a5c;
    }
    /* Result cards */
    .result-card {
        background: #F8FAFC;
        border: 1px solid #E2E8F0;
        border-radius: 8px;
        padding: 1.2rem;
        margin-bottom: 1rem;
    }
    /* DTSS score display */
    .dtss-score {
        font-size: 2.5rem;
        font-weight: 700;
        color: #0D2137;
    }
    /* Demo mode banner */
    .demo-banner {
        background: #EFF6FF;
        border-left: 4px solid #0D2137;
        padding: 0.8rem 1.2rem;
        border-radius: 0 6px 6px 0;
        font-size: 0.9rem;
        color: #0D2137;
        margin-bottom: 1rem;
    }
</style>
""", unsafe_allow_html=True)

st.markdown("""
<div class="saidock-header">
    <h1>SAIDock v1.0.0</h1>
    <p>Automated Drug-Target Docking and Assessment Pipeline | Pocket Detection | ADMET | ML Scoring | DTSS</p>
</div>
""", unsafe_allow_html=True)

if DEMO_MODE:
    st.markdown("""
<div class="demo-banner">
    <strong>Demo Mode</strong> — Pre-computed results are shown.
    Full docking requires a local AutoDock Vina installation.
    Clone from GitHub to run with your own HPC environment.
</div>
""", unsafe_allow_html=True)
# ─────────────────────────────────────────────────────────────────────

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
