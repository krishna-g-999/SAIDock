#!/usr/bin/env python3
"""
SAIDock Web Interface
Streamlit frontend for the SAIDock drug-target assessment pipeline.
"""
import sys, os
from pathlib import Path

#  Streamlit Cloud path fix 
_here = Path(__file__).parent.resolve()
for _p in [_here, _here.parent]:
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

#  Demo mode (True on Streamlit Cloud  no Vina binary) 
DEMO_MODE = (
    not Path('/home/krishnasalini-rs').exists() and
    not Path('/usr/local/bin/vina').exists() and
    not Path('/opt/vina').exists()
) or os.environ.get('SAIDOCK_DEMO', '').lower() == 'true'

# 
import streamlit as st
import requests
import json
import time

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

