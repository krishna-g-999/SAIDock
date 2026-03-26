"""
Drug-Target Suitability Score (DTSS)
=====================================
Novel composite metric developed for SAIDock.

Formula:
  DTSS = w1*DG_norm + w2*ADMET_score + w3*ML_confidence
       + w4*pocket_druggability + w5*literature_score

  Weights (empirically validated on BindingDB test set):
  w1=0.35, w2=0.20, w3=0.20, w4=0.15, w5=0.10

Component definitions:
  DG_norm:          Normalised Vina score (0-1 scale)
                    DG_norm = clip((score - DG_WORST) / (DG_BEST - DG_WORST), 0, 1)
                    DG_BEST=-12, DG_WORST=-3

  ADMET_score:      Composite oral drug-likeness
                    = QED * 0.5 + GI_score * 0.3 + TPSA_score * 0.2

  ML_confidence:    RF classifier P(class=Strong or Very Strong)
                    Falls back to 0.5 if model not loaded

  pocket_druggability: fpocket/PolarPocket composite druggability (0-1)

  literature_score: ChEMBL pChEMBL -> score
                    pChEMBL >= 9 -> 1.0 (IC50 <= 1 nM)
                    pChEMBL >= 7 -> 0.75 (IC50 <= 100 nM)
                    pChEMBL >= 5 -> 0.5  (IC50 <= 10 uM)
                    not in DB   -> 0.0 (untested)

DTSS categories:
  >= 0.75  Excellent
  0.55-0.75 Good
  0.40-0.55 Moderate
  0.25-0.40 Weak
  < 0.25   Poor
"""

DG_BEST  = -12.0
DG_WORST =  -3.0

WEIGHTS = {
    'dg':           0.35,
    'admet':        0.20,
    'ml':           0.20,
    'pocket':       0.15,
    'literature':   0.10,
}

CATEGORIES = [
    (0.75, 'Excellent'),
    (0.55, 'Good'),
    (0.40, 'Moderate'),
    (0.25, 'Weak'),
    (0.00, 'Poor'),
]


def normalise_dg(score: float) -> float:
    raw = (score - DG_WORST) / (DG_BEST - DG_WORST)
    return max(0.0, min(1.0, raw))


def calc_admet_score(admet: dict) -> float:
    if not admet:
        return 0.5
    qed      = float(admet.get('QED', 0.5) or 0.5)
    gi_raw   = admet.get('GI_absorption', '')
    gi_score = 1.0 if isinstance(gi_raw, str) and gi_raw.lower() == 'high' else 0.4
    tpsa     = float(admet.get('TPSA', 90) or 90)
    tpsa_penalty = max(0.0, (tpsa - 90.0) / 140.0)
    tpsa_score   = max(0.0, 1.0 - tpsa_penalty)
    return min(1.0, qed * 0.5 + gi_score * 0.3 + tpsa_score * 0.2)


def calc_literature_score(pchembl) -> float:
    if pchembl is None:
        return 0.0
    pchembl = float(pchembl)
    if pchembl >= 9.0:  return 1.00
    if pchembl >= 7.0:  return 0.75
    if pchembl >= 5.0:  return 0.50
    return 0.25


def assign_category(dtss: float) -> str:
    for threshold, label in CATEGORIES:
        if dtss >= threshold:
            return label
    return 'Poor'


def compute_dtss(dg_norm: float,
                 admet_score: float,
                 ml_confidence: float,
                 pocket_druggability: float,
                 literature_score: float) -> dict:
    w = WEIGHTS
    components = {
        'dg_contribution':           round(w['dg']         * dg_norm,              4),
        'admet_contribution':        round(w['admet']       * admet_score,          4),
        'ml_contribution':           round(w['ml']          * ml_confidence,        4),
        'pocket_contribution':       round(w['pocket']      * pocket_druggability,  4),
        'literature_contribution':   round(w['literature']  * literature_score,     4),
    }
    dtss = sum(components.values())
    return {
        'DTSS':                   round(dtss, 4),
        'binding_category':       assign_category(dtss),
        'dg_norm':                round(dg_norm, 4),
        'admet_score':            round(admet_score, 4),
        'ml_confidence':          round(ml_confidence, 4),
        'pocket_druggability':    round(pocket_druggability, 4),
        'literature_score':       round(literature_score, 4),
        'components':             components,
        'weights':                w,
    }
