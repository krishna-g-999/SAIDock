"""
SAIDock Report Generator — SAI-Net Design System
Inter font · Gold #F0A500 · Navy #040E1C gradient header
Professional text labels only — no emoji
"""
from pathlib import Path
from datetime import datetime

try:
    import plotly.graph_objects as go
    import plotly.io as pio
    HAS_PLOTLY = True
except ImportError:
    HAS_PLOTLY = False

# ── Runtime logo resolution (server-relative paths) ──────────────────
import glob as _glob

def _find_logo(pattern: str) -> str:
    """Return /static/<filename> if found, else empty string."""
    static = Path(__file__).parent.parent.parent / "web" / "static"
    hits = list(static.glob(pattern))
    return f"/static/{hits[0].name}" if hits else ""

# ── Design tokens ─────────────────────────────────────────────────────
GOLD       = "#F0A500"
GOLD_H     = "#D4920A"
NAVY1      = "#040E1C"
NAVY2      = "#112A47"
NAVY_MID   = "#0D2137"
NAVY_LT    = "#1A3A5C"
APP_BG     = "#EEF2F9"
CARD_BG    = "#FFFFFF"
CHART_BG   = "#FAFCFF"
TAB_BG     = "#F0F5FB"
BODY       = "#1E293B"
MUTED      = "#94A3B8"
BORDER     = "#CBD5E1"

# Evidence colours
E_HI_FG  = "#1B6B45"; E_HI_BG  = "#E6F5EE"
E_MO_FG  = "#9B5C00"; E_MO_BG  = "#FFF3E0"
E_LO_FG  = "#9B2335"; E_LO_BG  = "#FDE8EA"
E_ML_FG  = "#1D4ED8"; E_ML_BG  = "#DBEAFE"
E_ES_FG  = "#92710C"; E_ES_BG  = "#FFF3CD"

# Pocket chip colours
CHIP_MAP = {
    "Orthosteric": ("#1A3A5C", "#EEF4FF", "#C5D7F0"),
    "Allosteric":  ("#4A1F7A", "#F5F0FC", "#D8C8F5"),
    "Buried":      ("#0C4A6E", "#F0F9FF", "#BAE6FD"),
    "Cryptic":     ("#7A4000", "#FFF7EC", "#F2D9B0"),
    "Secondary":   ("#9B5C00", "#FFF3E0", "#F5D9A0"),
    "Surface":     ("#475569", "#F8FAFC", "#CBD5E1"),
}

CSS = """
@import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700;800&display=swap');

* { box-sizing: border-box; margin: 0; padding: 0; }

body {
  font-family: Inter, system-ui, -apple-system, sans-serif;
  background: #EEF2F9;
  color: #1E293B;
  font-size: 0.89rem;
  line-height: 1.65;
}

/* ── HEADER — exact SAIDock app match ── */
.site-header {
  background: linear-gradient(135deg, #040E1C 0%, #112A47 100%);
  border-bottom: 3px solid #F0A500;
  padding: 0 2.2rem;
  display: flex;
  align-items: center;
  justify-content: space-between;
  min-height: 70px;
  gap: 1.2rem;
}
.hdr-left {
  display: flex;
  align-items: center;
  gap: 1rem;
}
.hdr-logo-ring {
  width: 54px; height: 54px;
  border-radius: 50%;
  border: 2.5px solid #F0A500;
  box-shadow: 0 0 12px rgba(240,165,0,.40);
  background: #1A3A5C;
  display: flex; align-items: center; justify-content: center;
  font-size: 0.85rem; font-weight: 800; color: #F0A500;
  flex-shrink: 0;
  overflow: hidden;
}
.hdr-logo-ring img {
  width: 100%; height: 100%; object-fit: cover; border-radius: 50%;
}
.hdr-divider {
  width: 2px; height: 44px; flex-shrink: 0;
  background: linear-gradient(to bottom, transparent, #F0A500, transparent);
}
.hdr-brand { display: flex; flex-direction: column; }
.hdr-tool-name {
  font-size: 1.55rem;
  font-weight: 800;
  color: #F0A500;
  letter-spacing: -0.5px;
  line-height: 1.1;
}
.hdr-tool-name span {
  font-size: 0.82rem;
  font-weight: 500;
  color: #94A3B8;
  letter-spacing: 0;
}
.hdr-science {
  font-size: 0.72rem;
  font-weight: 700;
  color: #F0A500;
  letter-spacing: 1.4px;
  text-transform: uppercase;
  margin-top: 1px;
}
.hdr-sub {
  font-size: 0.68rem;
  font-weight: 400;
  color: #64748B;
  letter-spacing: 0.2px;
  margin-top: 2px;
}
.hdr-right {
  display: flex;
  align-items: center;
  gap: 14px;
}
.hdr-inst {
  text-align: right;
  line-height: 1.4;
}
.hdr-inst-name {
  font-size: 0.82rem;
  font-weight: 800;
  color: #CBD5E1;
  letter-spacing: 0.2px;
}
.hdr-inst-sub {
  font-size: 0.70rem;
  color: #64748B;
}
.hdr-inst-city {
  font-size: 0.68rem;
  font-weight: 700;
  color: #F0A500;
  letter-spacing: 0.8px;
  text-transform: uppercase;
  margin-top: 1px;
}
.hdr-sssihl-ring {
  width: 54px; height: 54px;
  border-radius: 50%;
  border: 2px solid rgba(255,255,255,.25);
  background: #fff;
  display: flex; align-items: center; justify-content: center;
  overflow: hidden;
  flex-shrink: 0;
}
.hdr-sssihl-ring img {
  width: 100%; height: 100%; object-fit: contain; padding: 3px;
}

/* ── Sub-nav ── */
.subnav {
  background: #0D2137;
  padding: 6px 2.2rem;
  font-size: 0.75rem;
  color: #64748B;
  border-bottom: 1px solid rgba(255,255,255,.06);
}
.subnav strong { color: #94A3B8; }

/* ── Page ── */
.page { max-width: 980px; margin: 0 auto; padding: 26px 20px 60px; }

/* ── Section title ── */
.stitle {
  font-size: 1.05rem;
  font-weight: 700;
  color: #0D2137;
  border-bottom: 2px solid #F0A500;
  padding-bottom: 6px;
  margin-bottom: 16px;
}

/* ── Cards ── */
.card {
  background: #fff;
  border: 1px solid #CBD5E1;
  border-radius: 10px;
  padding: 22px 24px;
  margin-bottom: 18px;
  box-shadow: 0 1px 4px rgba(4,14,28,.06);
}

/* ── DTSS hero ── */
.dtss-hero {
  display: flex; align-items: flex-start;
  gap: 28px; flex-wrap: wrap;
  padding: 20px; border-radius: 8px;
  border-left: 5px solid #F0A500;
  margin-bottom: 18px;
}
.dtss-num {
  font-size: 3.4rem; font-weight: 800;
  line-height: 1; letter-spacing: -1px;
}
.dtss-lbl {
  font-size: 0.68rem; font-weight: 800;
  text-transform: uppercase; letter-spacing: 1.6px;
  color: #94A3B8; margin-top: 4px;
}
.dtss-cat { font-size: 0.93rem; font-weight: 600; color: #0D2137; margin-top: 4px; }
.dtss-formula { font-size: 0.82rem; color: #475569; line-height: 1.9; }
.dtss-scale {
  font-size: 0.75rem; color: #64748B;
  margin-top: 8px; padding-top: 8px;
  border-top: 1px solid #E2E8F0;
}

/* ── Metric grid ── */
.mgrid {
  display: grid;
  grid-template-columns: repeat(auto-fit, minmax(128px, 1fr));
  gap: 12px; margin: 16px 0;
}
.mcard {
  background: #F0F5FB; border: 1px solid #CBD5E1;
  border-radius: 8px; padding: 14px 10px; text-align: center;
}
.mval { font-size: 1.25rem; font-weight: 800; color: #0D2137; }
.mlbl {
  font-size: 0.68rem; font-weight: 800;
  text-transform: uppercase; letter-spacing: 1.4px;
  color: #94A3B8; margin-top: 3px;
}

/* ── Tables ── */
table { width: 100%; border-collapse: collapse; font-size: 0.82rem; }
thead tr { background: linear-gradient(135deg, #040E1C, #112A47); }
thead th {
  padding: 10px 14px; text-align: left;
  color: #CBD5E1; font-size: 0.68rem; font-weight: 800;
  letter-spacing: 1px; text-transform: uppercase;
}
tbody tr:nth-child(even) { background: #F8FAFC; }
tbody tr:hover { background: #FFF8E8; transition: background .12s; }
td { padding: 9px 14px; border-bottom: 1px solid #F1F5F9; vertical-align: middle; }

/* ── Chips & badges ── */
.chip {
  display: inline-block; padding: 2px 9px;
  border-radius: 6px; font-size: 0.78rem; font-weight: 500;
  border-width: 1px; border-style: solid; white-space: nowrap;
}
.badge {
  display: inline-block; padding: 2px 10px;
  border-radius: 5px; font-size: 0.78rem; font-weight: 600;
}

/* ── Chart wrap ── */
.chart-wrap {
  margin: 14px 0; border: 1px solid #CBD5E1;
  border-radius: 8px; overflow: hidden; background: #FAFCFF;
}

/* ── About grid ── */
.agrid { display: grid; grid-template-columns: 1fr 1fr; gap: 16px; }
.ablock {
  background: #F0F5FB; border: 1px solid #CBD5E1;
  border-radius: 8px; padding: 16px 18px;
}
.ablock h4 {
  font-size: 0.68rem; font-weight: 800; color: #0D2137;
  text-transform: uppercase; letter-spacing: 1.6px;
  margin-bottom: 8px; padding-bottom: 5px;
  border-bottom: 1px solid #CBD5E1;
}
.ablock p, .ablock li { font-size: 0.82rem; color: #475569; line-height: 1.7; }
.ablock ul { padding-left: 14px; }

/* ── Footer ── */
.footer {
  text-align: center; font-size: 0.75rem; color: #94A3B8;
  margin-top: 36px; padding-top: 16px;
  border-top: 1px solid #CBD5E1;
}
.footer strong { color: #0D2137; }

@media print {
  body { background: #fff; }
  .site-header { -webkit-print-color-adjust: exact; print-color-adjust: exact; }
  .card { box-shadow: none; }
}
"""


def _chip(role: str) -> str:
    fg, bg, br = CHIP_MAP.get(role, ("#475569", "#F8FAFC", "#CBD5E1"))
    return (f'<span class="chip" style="color:{fg};background:{bg};'
            f'border-color:{br}">{role.upper()}</span>')


def _badge(text: str, level: str) -> str:
    lut = {"hi": (E_HI_FG, E_HI_BG), "mo": (E_MO_FG, E_MO_BG),
           "lo": (E_LO_FG, E_LO_BG), "ml": (E_ML_FG, E_ML_BG),
           "es": (E_ES_FG, E_ES_BG)}
    fg, bg = lut.get(level, ("#475569", "#F8FAFC"))
    return f'<span class="badge" style="color:{fg};background:{bg}">{text}</span>'


def _pf(ok) -> str:
    if ok is None:
        return f'<span style="color:{MUTED}">—</span>'
    return _badge("Pass", "hi") if ok else _badge("Fail", "lo")


def classify_pocket(pk: dict, dr_entry, all_dr: list) -> dict:
    comp   = pk.get("composite_score", pk.get("druggability_score", 0))
    vol    = pk.get("volume", 0)
    apol   = pk.get("apolar_sasa", 0)
    tot    = pk.get("total_sasa", 1)
    buried = apol / tot if tot > 0 else 0
    dg     = dr_entry.get("best_score", 0) if dr_entry else 0
    all_dg = [d.get("best_score", 0) for d in all_dr if d.get("best_score")]
    is_best = dg == min(all_dg) if all_dg else False

    if comp >= 0.70 and 400 <= vol <= 1400 and is_best:
        return dict(role="Orthosteric",
                    desc="Primary binding site — highest druggability and best docking affinity")
    if comp >= 0.65 and vol >= 700:
        return dict(role="Allosteric",
                    desc="High-volume secondary pocket — candidate for allosteric modulation")
    if buried >= 0.65 and comp >= 0.50:
        return dict(role="Buried",
                    desc="Deeply enclosed pocket — high apolar SASA fraction")
    if comp >= 0.50 and vol < 400:
        return dict(role="Cryptic",
                    desc="Small induced-fit pocket — may open upon ligand binding")
    if comp >= 0.40 and 400 <= vol <= 700:
        return dict(role="Secondary",
                    desc="Moderate druggability — candidate for fragment screening")
    return dict(role="Surface",
                desc="Low druggability — transient or solvent-exposed groove")


class ReportGenerator:
    def __init__(self, state: dict, outdir: str, log_fn=None):
        self.state  = state
        self.outdir = Path(outdir)
        self.outdir.mkdir(parents=True, exist_ok=True)
        self._log   = log_fn or (lambda msg, tag="": print(f"[{tag}] {msg}"))

    def log(self, msg, tag="INFO"):
        self._log(msg, tag)

    def generate(self):
        s      = self.state
        ml     = s.get("ml_results", {})
        adm    = s.get("admet", {})
        dr     = s.get("docking_results", [])
        pks    = s.get("pockets", [])
        ints   = s.get("interactions", {})
        target = s.get("target_id", "Unknown")
        ligand = s.get("ligand_name", "compound")
        nres   = s.get("n_residues", "—")
        date   = datetime.now().strftime("%Y-%m-%d %H:%M")

        dtss    = ml.get("DTSS", 0)
        cat     = ml.get("binding_category", "N/A")
        bestDG  = ml.get("best_docking_score", 0)
        admetSc = ml.get("admet_score", 0)
        mlConf  = ml.get("ml_confidence", 0.5)
        pktDrug = ml.get("pocket_druggability", 0)
        litSc   = ml.get("literature_score", 0)
        dgNorm  = ml.get("dg_norm", 0)
        mlAD    = "IN_DOMAIN" if ml.get("model_loaded") else "OUT_OF_DOMAIN"

        hb  = ints.get("hbonds", [])
        hy  = ints.get("hydrophobic", [])
        n_hb = len(hb) if isinstance(hb, list) else int(hb or 0)
        n_hy = len(hy) if isinstance(hy, list) else int(hy or 0)

        pk_map    = {p["pocket_id"]: p for p in pks}
        dr_sorted = sorted(dr, key=lambda x: x.get("best_score", 0))
        dr_map    = {r["pocket_id"]: r for r in dr_sorted}
        pk_class  = {pk["pocket_id"]: classify_pocket(pk, dr_map.get(pk["pocket_id"]), dr_sorted)
                     for pk in pks}

        # DTSS colour level
        if dtss >= 0.75:
            dc, db, dl = E_HI_FG,  E_HI_BG,  "hi"
        elif dtss >= 0.65:
            dc, db, dl = "#2D7A4F", "#F0FDF4", "hi"
        elif dtss >= 0.50:
            dc, db, dl = E_MO_FG,  E_MO_BG,  "mo"
        elif dtss >= 0.35:
            dc, db, dl = E_ES_FG,  E_ES_BG,  "es"
        else:
            dc, db, dl = E_LO_FG,  E_LO_BG,  "lo"

        ml_badge = (_badge("IN DOMAIN", "hi") if mlAD == "IN_DOMAIN"
                    else _badge("OUT OF DOMAIN", "es"))

        # Charts
        dtss_chart = self._dtss_chart(ml, dtss) if HAS_PLOTLY else ""
        vina_chart = self._vina_chart(dr_sorted, pk_class) if HAS_PLOTLY else ""
        drug_chart = self._drug_chart(pks, pk_class) if HAS_PLOTLY else ""

        # Docking rows
        dock_rows = ""
        for i, r in enumerate(dr_sorted):
            pid   = r.get("pocket_id")
            dg_v  = r.get("best_score", 0)
            m3    = r.get("mean_top3", 0)
            poses = r.get("n_poses", 9)
            eng   = r.get("engine", "AutoDock Vina")
            drug  = r.get("pocket_drugg_score",
                    pk_map.get(pid, {}).get("composite_score",
                    pk_map.get(pid, {}).get("druggability_score", 0)))
            role  = pk_class.get(pid, {}).get("role", "Surface")
            best  = i == 0
            dg_c  = E_HI_FG if dg_v <= -8 else E_MO_FG if dg_v <= -6 else E_LO_FG
            lbl   = "Best Pocket &mdash; " if best else ""
            dock_rows += f"""
      <tr style="background:{'#FFFBEB' if best else 'transparent'}">
        <td style="font-weight:{'700' if best else '500'};color:{'#92400E' if best else BODY}">{lbl}Pocket {pid}</td>
        <td style="font-weight:700;color:{dg_c}">{dg_v:.3f}</td>
        <td style="color:#475569">{m3:.3f}</td>
        <td>{poses}</td>
        <td style="color:#475569">{eng}</td>
        <td style="font-weight:600;color:{'#1B6B45' if drug>=0.7 else E_MO_FG if drug>=0.4 else MUTED}">{drug:.3f}</td>
        <td>{_chip(role)}</td>
      </tr>"""

        # Pocket explorer rows
        pocket_rows = ""
        for pk in sorted(pks, key=lambda x: -x.get("composite_score", 0)):
            pid   = pk.get("pocket_id")
            vol   = pk.get("volume", 0)
            sasa  = pk.get("total_sasa", 0)
            hydro = pk.get("hydrophobicity", 0)
            drug2 = pk.get("composite_score", pk.get("druggability_score", 0))
            nsph  = pk.get("n_spheres", 0)
            apol  = pk.get("apolar_sasa", 0)
            bur   = (apol / sasa * 100) if sasa > 0 else 0
            dr_e  = dr_map.get(pid)
            dg2   = dr_e.get("best_score", 0) if dr_e else 0
            cls   = pk_class.get(pid, {})
            role  = cls.get("role", "Surface")
            desc  = cls.get("desc", "")
            dg_c  = E_HI_FG if dg2 <= -8 else E_MO_FG if dg2 <= -6 else E_LO_FG
            pocket_rows += f"""
      <tr title="{desc}">
        <td style="font-weight:600">Pocket {pid}</td>
        <td style="font-weight:700;color:{dg_c}">{dg2:.3f}</td>
        <td>{vol:.1f}</td><td>{sasa:.1f}</td>
        <td>{bur:.1f}%</td><td>{hydro:.2f}</td><td>{int(nsph)}</td>
        <td style="font-weight:600;color:{'#1B6B45' if drug2>=0.7 else E_MO_FG if drug2>=0.4 else MUTED}">{drug2:.3f}</td>
        <td>{_chip(role)}</td>
      </tr>"""

        # ADMET rows
        PROPS = [
            ("MW",                 "Molecular Weight",              "g/mol",  "&le; 500",         lambda v: float(v)<=500),
            ("LogP",               "Lipophilicity",                 "LogP",   "&le; 5",            lambda v: float(v)<=5),
            ("HBD",                "H-Bond Donors",                 "",       "&le; 5",            lambda v: float(v)<=5),
            ("HBA",                "H-Bond Acceptors",              "",       "&le; 10",           lambda v: float(v)<=10),
            ("TPSA",               "Topological PSA",               "A²",     "&le; 140",          lambda v: float(v)<=140),
            ("RotBonds",           "Rotatable Bonds",               "",       "&le; 10",           lambda v: float(v)<=10),
            ("HeavyAtoms",         "Heavy Atom Count",              "",       "20&ndash;70",       lambda v: 20<=float(v)<=70),
            ("MolRefractivity",    "Molar Refractivity",            "",       "40&ndash;130",      lambda v: 40<=float(v)<=130),
            ("AromaticRings",      "Aromatic Ring Count",           "",       "&le; 4",            lambda v: float(v)<=4),
            ("QED",                "Drug-likeness QED",             "",       "&ge; 0.5",          lambda v: float(v)>=0.5),
            ("FractionCSP3",       "Fsp3",                          "",       "&ge; 0.25",         lambda v: float(v)>=0.25),
            ("GI_absorption",      "GI Absorption",                 "",       "High",              lambda v: "high" in str(v).lower()),
            ("BBB_penetrant",      "BBB Penetrant",                 "",       "—",                 None),
            ("Pgp_substrate",      "P-gp Substrate",                "",       "Unlikely",          lambda v: "unlikely" in str(v).lower()),
            ("Solubility_class",   "Aqueous Solubility",            "",       "—",                 None),
            ("LogSw_estimated",    "LogSw ESOL",                    "",       "&gt; &minus;4",     lambda v: float(v)>-4),
            ("Lipinski_violations","Lipinski Violations",           "",       "&le; 1",            lambda v: float(v)<=1),
            ("Veber_pass",         "Veber Oral Bioavailability",    "",       "Pass",              lambda v: str(v).lower() in ["true","1","pass"]),
            ("Ghose_pass",         "Ghose Filter",                  "",       "Pass",              lambda v: str(v).lower() in ["true","1","pass"]),
        ]
        admet_rows = ""
        for key, label, unit, thresh, test in PROPS:
            val = adm.get(key)
            if val is None: continue
            status = _pf(None) if test is None else _pf(test(val) if True else None)
            try:
                if test is not None: status = _pf(test(val))
            except: pass
            disp = f"{float(val):.3f}" if isinstance(val, (int, float)) else str(val)
            if unit: disp += f' <span style="color:{MUTED};font-size:0.75rem">{unit}</span>'
            admet_rows += f"""
      <tr>
        <td style="font-weight:500">{label}</td>
        <td style="font-weight:600">{disp}</td>
        <td style="color:{MUTED}">{thresh}</td>
        <td>{status}</td>
      </tr>"""

        ml_rows_html = "".join(
            f'<tr><td style="font-weight:500">{k}</td><td style="font-weight:600">{v}</td></tr>'
            for k, v in [
                ("DTSS Score",               f"{dtss:.4f}"),
                ("Binding Category",         cat),
                ("DG Normalised",            f"{dgNorm:.4f}"),
                ("ADMET Score",              f"{admetSc:.4f}"),
                ("ML Confidence",            f"{mlConf:.4f}"),
                ("Pocket Druggability",      f"{pktDrug:.4f}"),
                ("Literature Prior Score",   f"{litSc:.4f}"),
                ("ML Applicability Domain",  mlAD),
                ("ChEMBL pChEMBL",           ml.get("chembl_pchembl","Not in database")),
            ])

        hbond_rows = ""
        if isinstance(hb, list):
            for bond in hb[:25]:
                if isinstance(bond, dict):
                    aa = bond.get("residue_name",""); res = bond.get("residue_id","")
                    ch = bond.get("chain","A"); d = bond.get("distance", bond.get("dist",""))
                    t  = bond.get("type","Hydrogen Bond")
                    dd = f"{d:.2f} A" if isinstance(d, float) else str(d)
                    hbond_rows += (f"<tr><td>{aa} {res} (Chain {ch})</td>"
                                   f"<td>{t}</td><td>{dd}</td></tr>")

        inter_card = ""
        if hbond_rows:
            inter_card = f"""
<div class="card">
  <div class="stitle">Key Binding Interactions &mdash; Best Pocket</div>
  <p style="font-size:0.82rem;color:#475569;margin-bottom:12px">
    Hydrogen bonds: <strong>{n_hb}</strong> &nbsp;|&nbsp;
    Hydrophobic contacts: <strong>{n_hy}</strong>
  </p>
  <table>
    <thead><tr><th>Residue</th><th>Interaction Type</th><th>Distance</th></tr></thead>
    <tbody>{hbond_rows}</tbody>
  </table>
</div>"""

        # Resolve logos at render time
        sai_src   = _find_logo("[Ss][Aa][Ii]*[Nn][Ee][Tt]*.png") or _find_logo("sai*.png")
        ssihl_src = _find_logo("[Ss][Ss][Ss][Ii][Hh][Ll]*.png") or _find_logo("sssihl*.png")

        if sai_src:
            left_logo = f'<div class="hdr-logo-ring"><img src="{sai_src}" alt="SAI-Net"></div>'
        else:
            left_logo = '<div class="hdr-logo-ring">SD</div>'

        if ssihl_src:
            right_logo = f'<div class="hdr-sssihl-ring"><img src="{ssihl_src}" alt="SSSIHL"></div>'
        else:
            right_logo = '<div class="hdr-sssihl-ring" style="font-size:0.52rem;font-weight:800;color:#1B2A4A;text-align:center;padding:4px;line-height:1.2">SSSIHL</div>'

        plotly_cdn = '<script src="https://cdn.plot.ly/plotly-2.26.0.min.js"></script>'

        html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>SAIDock Report &mdash; {ligand} vs {target}</title>
{plotly_cdn if HAS_PLOTLY else ''}
<style>{CSS}</style>
</head>
<body>

<header class="site-header">
  <div class="hdr-left">
    {left_logo}
    <div class="hdr-divider"></div>
    <div class="hdr-brand">
      <div class="hdr-tool-name">SAIDock <span>V 1.0</span></div>
      <div class="hdr-science">Science for Society</div>
      <div class="hdr-sub">
        Automated Drug-Target Docking &amp; Assessment
        &nbsp;&middot;&nbsp; NDD Research Platform
        &nbsp;&middot;&nbsp; SAI-Net Framework
      </div>
    </div>
  </div>
  <div class="hdr-right">
    <div class="hdr-inst">
      <div class="hdr-inst-name">SSSIHL</div>
      <div class="hdr-inst-sub">Sri Sathya Sai Institute<br>of Higher Learning</div>
      <div class="hdr-inst-city">Puttaparthi &nbsp;&middot;&nbsp; Est. 1981</div>
    </div>
    {right_logo}
  </div>
</header>

<div class="subnav">
  <strong>Ligand:</strong> {ligand}
  &nbsp;&nbsp;<strong>Target:</strong> {target} ({nres} residues)
  &nbsp;&nbsp;<strong>Generated:</strong> {date}
  &nbsp;&nbsp;SAIDock V 1.0 &nbsp;&middot;&nbsp; SAI-Net Research Framework
</div>

<div class="page">

<div class="card">
  <div class="dtss-hero" style="background:{db}">
    <div>
      <div class="dtss-num" style="color:{dc}">{dtss:.3f}</div>
      <div class="dtss-lbl">Drug-Target Suitability Score</div>
      <div class="dtss-cat">{cat}</div>
      <div style="margin-top:10px">{ml_badge}</div>
    </div>
    <div style="flex:1;min-width:260px">
      <div class="stitle" style="font-size:0.93rem;margin-bottom:10px">DTSS Scoring Formula</div>
      <div class="dtss-formula">
        DTSS = 0.35 &times; &Delta;G<sub>norm</sub>
        + 0.20 &times; ADMET
        + 0.20 &times; ML<sub>confidence</sub>
        + 0.15 &times; Pocket<sub>druggability</sub>
        + 0.10 &times; Literature<sub>prior</sub>
      </div>
      <div class="dtss-scale">
        &ge; 0.75 Strong &nbsp;|&nbsp;
        0.65&ndash;0.75 Good &nbsp;|&nbsp;
        0.50&ndash;0.65 Moderate &nbsp;|&nbsp;
        0.35&ndash;0.50 Weak &nbsp;|&nbsp;
        &lt; 0.35 Poor
      </div>
    </div>
  </div>
  <div class="mgrid">
    <div class="mcard"><div class="mval" style="color:{GOLD}">{bestDG:.3f}</div><div class="mlbl">Best &Delta;G kcal/mol</div></div>
    <div class="mcard"><div class="mval">{len(dr)}</div><div class="mlbl">Pockets Screened</div></div>
    <div class="mcard"><div class="mval">{n_hb}</div><div class="mlbl">Hydrogen Bonds</div></div>
    <div class="mcard"><div class="mval">{n_hy}</div><div class="mlbl">Hydrophobic</div></div>
    <div class="mcard"><div class="mval">{admetSc:.3f}</div><div class="mlbl">ADMET Score</div></div>
    <div class="mcard"><div class="mval">{pktDrug:.3f}</div><div class="mlbl">Pocket Druggability</div></div>
  </div>
  {f'<div class="chart-wrap">{dtss_chart}</div>' if dtss_chart else ''}
</div>

<div class="card">
  <div class="stitle">Multi-Pocket Docking Results</div>
  {f'<div class="chart-wrap" style="margin-bottom:16px">{vina_chart}</div>' if vina_chart else ''}
  <table>
    <thead><tr>
      <th>Binding Pocket</th><th>Best &Delta;G (kcal/mol)</th>
      <th>Mean Top-3</th><th>Poses</th><th>Engine</th>
      <th>Druggability</th><th>Site Classification</th>
    </tr></thead>
    <tbody>{dock_rows}</tbody>
  </table>
</div>

<div class="card">
  <div class="stitle">Binding Pocket Explorer</div>
  <p style="font-size:0.82rem;color:#475569;margin-bottom:14px">
    All fpocket-detected binding sites ranked by composite druggability score.
    Hover a chip for scientific rationale.
  </p>
  {f'<div class="chart-wrap" style="margin-bottom:16px">{drug_chart}</div>' if drug_chart else ''}
  <table>
    <thead><tr>
      <th>Pocket</th><th>Best &Delta;G</th><th>Volume (A³)</th>
      <th>SASA (A²)</th><th>Buried %</th><th>Hydrophobicity</th>
      <th>Alpha Spheres</th><th>Druggability</th><th>Classification</th>
    </tr></thead>
    <tbody>{pocket_rows}</tbody>
  </table>
</div>

{inter_card}

<div class="card">
  <div class="stitle">ADMET Drug-Likeness Assessment</div>
  <table>
    <thead><tr><th>Property</th><th>Value</th><th>Threshold</th><th>Status</th></tr></thead>
    <tbody>{admet_rows}</tbody>
  </table>
</div>

<div class="card">
  <div class="stitle">ML Scoring Components</div>
  <table>
    <thead><tr><th>Component</th><th>Value</th></tr></thead>
    <tbody>{ml_rows_html}</tbody>
  </table>
</div>

<div class="card">
  <div class="stitle">Computational Methods</div>
  <table>
    <thead><tr><th>Stage</th><th>Method / Software</th><th>Parameters</th></tr></thead>
    <tbody>
      <tr><td style="font-weight:500">Target Preparation</td><td>Biopython 1.81 + pdbfixer + OpenBabel 3.1</td><td>pH 7.4; chain A; Gasteiger charges</td></tr>
      <tr><td style="font-weight:500">Pocket Detection</td><td>fpocket 4.0 (Le Guilloux et al., 2009)</td><td>Alpha sphere 3.0&ndash;7.5 A; composite druggability score</td></tr>
      <tr><td style="font-weight:500">Ligand Preparation</td><td>RDKit ETKDG v3 + MMFF94</td><td>3D conformer; torsion minimisation; Gasteiger charges</td></tr>
      <tr><td style="font-weight:500">Molecular Docking</td><td>AutoDock Vina 1.2 (Eberhardt et al., 2021)</td><td>Exhaustiveness = 32; poses = 9; multi-pocket mode</td></tr>
      <tr><td style="font-weight:500">Interaction Analysis</td><td>ProLIF + geometry-based cutoffs</td><td>H-bond &le; 3.5 A; hydrophobic &le; 4.5 A</td></tr>
      <tr><td style="font-weight:500">ADMET Profiling</td><td>RDKit 2023.09 descriptors</td><td>Lipinski Ro5; Veber; Ghose; BBB; GI absorption; ESOL</td></tr>
      <tr><td style="font-weight:500">ML Scoring</td><td>Random Forest + DTSS composite</td><td>Morgan FP-512 + 12 descriptors; LOTO-CV; applicability domain</td></tr>
      <tr><td style="font-weight:500">Literature Prior</td><td>ChEMBL 33 REST API</td><td>IC&#8325;&#8320; / Kd / Ki &rarr; pChEMBL &rarr; literature_score</td></tr>
      <tr><td style="font-weight:500">Pocket Classification</td><td>Schmidtke &amp; Barril, JCTC 2010</td><td>Volume; apolar SASA; composite druggability; docking rank</td></tr>
    </tbody>
  </table>
</div>

<div class="card">
  <div class="stitle">About SAIDock</div>
  <div class="agrid">
    <div class="ablock">
      <h4>Overview</h4>
      <p>SAIDock is an automated multi-stage pipeline for structure-based
      drug-target interaction assessment. It integrates pocket detection,
      molecular docking, ADMET profiling, ML scoring, and literature
      validation into a single reproducible workflow. The DTSS ranks
      compound-target pairs beyond raw docking affinity alone.</p>
    </div>
    <div class="ablock">
      <h4>Key Capabilities</h4>
      <ul>
        <li>Automated multi-pocket detection and druggability classification</li>
        <li>AutoDock Vina docking across all detected pockets</li>
        <li>19-parameter ADMET profiling via RDKit descriptors</li>
        <li>Random Forest ML scoring with applicability domain</li>
        <li>Real-time ChEMBL literature prior integration</li>
        <li>DTSS composite scoring for lead prioritisation</li>
      </ul>
    </div>
    <div class="ablock">
      <h4>Developed By</h4>
      <p>Department of Bioinformatics<br>
      Sri Sathya Sai Institute of Higher Learning (SSSIHL)<br>
      Puttaparthi, Andhra Pradesh, India &mdash; Est. 1981<br><br>
      Part of the SAI-Net Research Framework for neurodegenerative
      disease computational research.</p>
    </div>
    <div class="ablock">
      <h4>Citation</h4>
      <p>For academic research use only. Not intended for clinical
      diagnosis or therapeutic decision-making.<br><br>
      If using SAIDock in published work, please cite:<br>
      <em>Gunanathan K. et al., SAIDock: Automated Drug-Target
      Suitability Assessment (2026), SSSIHL, Puttaparthi.</em></p>
    </div>
  </div>
</div>

<div class="footer">
  <strong>SAIDock V 1.0</strong>
  &nbsp;&middot;&nbsp; Sri Sathya Sai Institute of Higher Learning, Puttaparthi
  &nbsp;&middot;&nbsp; {date}
  &nbsp;&middot;&nbsp; SAI-Net Research Framework
  &nbsp;&middot;&nbsp; <em>For academic research use only</em>
</div>

</div>
</body>
</html>"""

        out = self.outdir / "saidock_report.html"
        out.write_text(html, encoding="utf-8")
        self.log(f"HTML report saved: {out}", "OK")
        return str(out)

    def _cl(self, title, **kw):
        return dict(
            title=dict(text=title, font=dict(size=13, family="Inter,sans-serif", color=NAVY_MID)),
            plot_bgcolor=CHART_BG, paper_bgcolor=CHART_BG,
            font=dict(family="Inter,sans-serif", size=11, color=BODY), **kw)

    def _dtss_chart(self, ml, dtss):
        try:
            W = {"DG (x0.35)": 0.35*ml.get("dg_norm",0),
                 "ADMET (x0.20)": 0.20*ml.get("admet_score",0),
                 "ML conf (x0.20)": 0.20*ml.get("ml_confidence",0.5),
                 "Pocket (x0.15)": 0.15*ml.get("pocket_druggability",0),
                 "Literature (x0.10)": 0.10*ml.get("literature_score",0)}
            fig = go.Figure(go.Bar(
                x=list(W.values()), y=list(W.keys()), orientation="h",
                marker_color=[NAVY_LT, E_HI_FG, "#7C3AED", GOLD, MUTED],
                text=[f"{v:.3f}" for v in W.values()], textposition="outside",
                hovertemplate="%{y}: %{x:.4f}<extra></extra>"))
            fig.add_vline(x=dtss, line_dash="dash", line_color=E_LO_FG, line_width=1.5,
                          annotation_text=f"DTSS={dtss:.3f}",
                          annotation_font_color=E_LO_FG, annotation_font_size=10,
                          annotation_position="top right")
            fig.update_layout(**self._cl("DTSS Component Breakdown",
                xaxis_title="Weighted contribution",
                xaxis=dict(range=[0, max(dtss*1.18,0.4)], gridcolor="#E2E8F0"),
                margin=dict(l=165,r=70,t=46,b=38), height=268))
            return pio.to_html(fig, include_plotlyjs=False, full_html=False,
                               config={"displayModeBar": False})
        except Exception as e:
            return f"<!-- dtss chart error: {e} -->"

    def _vina_chart(self, dr_sorted, pk_class):
        try:
            CMAP = {"Orthosteric":NAVY_LT,"Allosteric":"#4A1F7A",
                    "Buried":"#0C4A6E","Cryptic":"#7A4000",
                    "Secondary":E_MO_FG,"Surface":"#475569"}
            names  = [f"Pocket {r['pocket_id']}" for r in dr_sorted]
            vals   = [r.get("best_score",0) for r in dr_sorted]
            errs   = [abs(r.get("best_score",0)-r.get("mean_top3",0)) for r in dr_sorted]
            cols   = [CMAP.get(pk_class.get(r["pocket_id"],{}).get("role","Surface"), NAVY_MID)
                      for r in dr_sorted]
            fig = go.Figure(go.Bar(
                x=vals, y=names, orientation="h", marker_color=cols,
                error_x=dict(type="data",array=errs,color="#94A3B8",thickness=1.5,width=4),
                text=[f"{v:.2f}" for v in vals], textposition="outside",
                hovertemplate="Pocket %{y}: %{x:.3f} kcal/mol<extra></extra>"))
            fig.update_layout(**self._cl("Binding Affinity by Pocket — AutoDock Vina",
                xaxis_title="Binding Affinity (kcal/mol)",
                xaxis=dict(gridcolor="#E2E8F0"),
                margin=dict(l=110,r=80,t=46,b=38),
                height=max(220, 50+50*len(dr_sorted))))
            return pio.to_html(fig, include_plotlyjs=False, full_html=False,
                               config={"displayModeBar": False})
        except Exception as e:
            return f"<!-- vina chart error: {e} -->"

    def _drug_chart(self, pks, pk_class):
        try:
            ps = sorted(pks, key=lambda x: -x.get("composite_score",0))
            names  = [f"Pocket {p['pocket_id']}" for p in ps]
            drug_s = [p.get("composite_score", p.get("druggability_score",0)) for p in ps]
            vols   = [min(p.get("volume",0)/1500,1.0) for p in ps]
            hydros = [min(p.get("hydrophobicity",0)/80,1.0) for p in ps]
            fig = go.Figure()
            fig.add_trace(go.Bar(name="Druggability Score", x=names, y=drug_s,
                                 marker_color=GOLD, text=[f"{v:.3f}" for v in drug_s],
                                 textposition="outside"))
            fig.add_trace(go.Bar(name="Volume / 1500", x=names, y=vols,
                                 marker_color=NAVY_MID, opacity=0.60))
            fig.add_trace(go.Bar(name="Hydrophobicity / 80", x=names, y=hydros,
                                 marker_color=E_HI_FG, opacity=0.65))
            fig.add_hline(y=0.7, line_dash="dot", line_color=E_LO_FG, line_width=1.2,
                          annotation_text="Druggable threshold 0.70",
                          annotation_font_color=E_LO_FG, annotation_font_size=10)
            fig.update_layout(**self._cl("Pocket Druggability Profile",
                yaxis_title="Normalised score",
                yaxis=dict(gridcolor="#E2E8F0"), barmode="group",
                legend=dict(orientation="h", y=1.12, x=0.5, xanchor="center",
                            font=dict(size=10)),
                margin=dict(l=55,r=55,t=60,b=38), height=305))
            return pio.to_html(fig, include_plotlyjs=False, full_html=False,
                               config={"displayModeBar": False})
        except Exception as e:
            return f"<!-- drug chart error: {e} -->"
