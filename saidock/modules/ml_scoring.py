#!/usr/bin/env python3
"""
Stage 07: ML Scoring + DTSS
- Trains per-target Random Forest classifiers on ChEMBL activity data
- Predicts binding probability (ML confidence) for query compound
- Computes Drug-Target Suitability Score (DTSS)
"""
import json
import pickle
import numpy as np
from pathlib import Path

# RDKit
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors

# scikit-learn
from sklearn.ensemble import RandomForestClassifier, GradientBoostingClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.metrics import roc_auc_score
import warnings
warnings.filterwarnings('ignore')

# ── DTSS constants ──────────────────────────────────────────────────────────
# Standard weights (data-rich targets: >=100 training samples)
DTSS_WEIGHTS = dict(dg=0.35, admet=0.20, ml=0.20, pocket=0.15, lit=0.10)

# Reduced ML weight for data-scarce targets (<50 training samples).
# Docking score gets higher weight as ML evidence is limited.
DTSS_WEIGHTS_SPARSE = dict(dg=0.45, admet=0.20, ml=0.10, pocket=0.15, lit=0.10)

# RNA-binding protein targets — RNA displacement score included in DTSS
RNA_BINDING_TARGETS = {"TDP43_RRM1", "TDP43_RRM2", "TDP43_pept"}
RNA_SITE_TARGETS    = {"5HDN", "6N37", "2N3X", "2N3Y", "4IUF",
                        "4BS2", "2KFN", "6CFH", "6CZ8", "2C9V"}


# Targets trained on <50 samples — use sparse weights
# CDK6 and CK2a added: natural product library is outside
# ATP-competitive inhibitor training domain (applicability domain limit)
SPARSE_TARGETS = {'TDP43_RRM1', 'TDP43_pept', 'Parkin_UBL', 'HTT', 'CDK6', 'CK2a'  # AD mismatch with natural product library
}
DG_BEST      = -12.0    # kcal/mol  (strongest realistic binder)
DG_WORST     =  -3.0    # kcal/mol  (non-binder threshold)

# Activity thresholds (pChEMBL)
ACTIVE_THRESHOLD   = 6.0   # pChEMBL >= 6  → active (IC50 <= 1 µM)
INACTIVE_THRESHOLD = 5.0   # pChEMBL <  5  → inactive

FP_BITS      = 2048
FP_RADIUS    = 2


class MLScorer:
    """
    Per-target Random Forest ML scorer.

    Usage:
        # Train
        scorer = MLScorer(model_dir='models/')
        scorer.train(chembl_activities, target_id='CK2a')

        # Predict
        scorer = MLScorer(model_dir='models/')
        result = scorer.score(smiles, docking_results, admet, pockets, lit_score)
    """

    def __init__(self, model_dir: str = 'models/', logger=None,
                 ligand_smiles: str = None, target_id: str = None,
                 docking_results: list = None, admet: dict = None,
                 pockets: list = None, lit_score: float = None,
                 outdir: str = None, **kwargs):
        self.model_dir       = Path(model_dir)
        self.model_dir.mkdir(parents=True, exist_ok=True)
        self.logger          = logger
        self._rf             = None
        self._scaler         = None
        self._target_id      = target_id
        self._ligand_smiles  = ligand_smiles
        self._docking_results= docking_results or []
        self._admet          = admet or {}
        self._pockets        = pockets or []
        self._lit_score      = lit_score if lit_score is not None else 0.0


    def log(self, msg, level='INFO'):
        if self.logger:
            self.logger.log(msg, level)

    # ── Training ──────────────────────────────────────────────────────────────
    def train(self, activities: list, target_id: str,
              n_estimators: int = 300,
              min_samples: int = 5) -> dict:
        """
        Train a RF classifier on ChEMBL activity data.

        activities : list of dicts with keys:
            smiles, pchembl_value, dg_estimated
        target_id  : string label used to save the model
        Returns    : training metrics dict
        """
        self._target_id = target_id
        self.log(f'Training ML model for {target_id} | {len(activities)} records')

        X, y, valid_smiles = self._featurise_activities(activities)
        if len(X) < min_samples:
            self.log(
                f'  Only {len(X)} valid molecules — need >= {min_samples}. '
                f'Skipping training.', 'WARN'
            )
            return {'status': 'skipped', 'n_samples': len(X)}

        class_counts = np.bincount(y)
        self.log(
            f'  Active: {class_counts[1] if len(class_counts)>1 else 0}  '
            f'Inactive: {class_counts[0]}  Total: {len(y)}'
        )

        # Scale features
        self._scaler = StandardScaler()
        Xs = self._scaler.fit_transform(X)

        # Train RF
        self._rf = RandomForestClassifier(
            n_estimators  = n_estimators,
            max_depth     = 10,
            min_samples_leaf = 3,
            class_weight  = 'balanced',
            n_jobs        = -1,
            random_state  = 42,
        )

        # ── Cross-validation with multiple metrics ──────────────────────────
        from sklearn.metrics import (average_precision_score,
                                     matthews_corrcoef, roc_auc_score,
                                     precision_recall_curve)
        from sklearn.model_selection import StratifiedKFold

        n_splits   = max(2, min(5, int(min(class_counts)))) if len(class_counts)>1 and min(class_counts)>=2 else 2
        cv         = StratifiedKFold(n_splits=n_splits, shuffle=True,
                                     random_state=42)
        roc_aucs, prc_aucs, mccs = [], [], []
        try:
            for train_idx, val_idx in cv.split(Xs, y):
                Xtr, Xval = Xs[train_idx], Xs[val_idx]
                ytr, yval = y[train_idx],  y[val_idx]
                self._rf.fit(Xtr, ytr)
                proba = self._rf.predict_proba(Xval)[:, 1]
                pred  = self._rf.predict(Xval)
                roc_aucs.append(roc_auc_score(yval, proba))
                prc_aucs.append(average_precision_score(yval, proba))
                mccs.append(matthews_corrcoef(yval, pred))
        except Exception as e:
            self.log(f'  Cross-validation failed: {e}', 'WARN')
            roc_aucs = [0.0]; prc_aucs = [0.0]; mccs = [0.0]

        # ── Scaffold split evaluation ─────────────────────────────────────────
        scaffold_auc = self._scaffold_split_eval(X, y, valid_smiles)

        # ── Fit final model on ALL data ───────────────────────────────────────
        self._rf.fit(Xs, y)

        # Feature importance
        desc_names  = self._descriptor_names()
        importances = self._rf.feature_importances_
        top5_idx    = np.argsort(importances)[::-1][:5]
        top5 = [(desc_names[i] if i < len(desc_names) else f'FP_{i}',
                 round(float(importances[i]), 4))
                for i in top5_idx]

        # Class balance warning
        imbalance = max(class_counts) / max(min(class_counts), 1)
        if imbalance > 10:
            self.log(
                f'  WARNING: class imbalance = {imbalance:.1f}:1 — '
                f'ROC-AUC is optimistic. Report AUPRC as primary metric.', 'WARN'
            )

        metrics = {
            'target_id':         target_id,
            'n_samples':         len(y),
            'n_active':          int(class_counts[1]) if len(class_counts)>1 else 0,
            'n_inactive':        int(class_counts[0]),
            'class_imbalance':   round(float(imbalance), 2),
            'cv_roc_auc_mean':   round(float(np.mean(roc_aucs)), 4),
            'cv_roc_auc_std':    round(float(np.std(roc_aucs)),  4),
            'cv_auprc_mean':     round(float(np.mean(prc_aucs)), 4),
            'cv_auprc_std':      round(float(np.std(prc_aucs)),  4),
            'cv_mcc_mean':       round(float(np.mean(mccs)),     4),
            'cv_mcc_std':        round(float(np.std(mccs)),      4),
            'scaffold_roc_auc':  round(float(scaffold_auc),      4),
            'auprc_baseline':    round(float(min(class_counts)/len(y)), 4),
            'top5_features':     top5,
            'status':            'trained',
            'note': (
                'ROC-AUC is reported for comparison but AUPRC is the primary metric '
                'due to class imbalance. Scaffold-split AUC tests generalisation '
                'to structurally dissimilar compounds.'
            ),
        }

        self.log(
            f'  ROC-AUC = {np.mean(roc_aucs):.3f} ± {np.std(roc_aucs):.3f} | '
            f'AUPRC = {np.mean(prc_aucs):.3f} (baseline={min(class_counts)/len(y):.3f}) | '
            f'MCC = {np.mean(mccs):.3f}', 'OK'
        )
        self.log(
            f'  Scaffold-split AUC = {scaffold_auc:.3f} | '
            f'Top feature: {top5[0][0]} ({top5[0][1]})'
        )

        # Save
        model_path   = self.model_dir / f'{target_id}_rf.pkl'
        scaler_path  = self.model_dir / f'{target_id}_scaler.pkl'
        metrics_path = self.model_dir / f'{target_id}_metrics.json'

        with open(model_path,   'wb') as f: pickle.dump(self._rf,      f)
        with open(scaler_path,  'wb') as f: pickle.dump(self._scaler,  f)
        with open(metrics_path, 'w')  as f: json.dump(metrics,          f, indent=2)

        self.log(f'  Model saved: {model_path}', 'OK')
        return metrics

    # ── Prediction + DTSS ─────────────────────────────────────────────────────
    def score(self, smiles: str = None,
              docking_results: list = None,
              admet: dict = None,
              pockets: list = None,
              lit_score: float = None,
              target_id: str = None) -> dict:
        """
        Compute ML confidence + DTSS for a query compound.
        All parameters optional — falls back to values stored in __init__
        when called from pipeline as ml.score() with no arguments.
        """
        # Resolve from stored state if not passed directly
        smiles          = smiles          or self._ligand_smiles  or ''
        docking_results = docking_results or self._docking_results or []
        admet           = admet           or self._admet          or {}
        pockets         = pockets         or self._pockets        or []
        target_id       = target_id       or self._target_id
        if lit_score is None:
            lit_score   = self._lit_score if self._lit_score is not None else 0.0
        # Load model if available
        ml_confidence = 0.50   # neutral default (model not trained yet)
        model_loaded  = False

        if target_id:
            ml_confidence, model_loaded = self._predict(smiles, target_id)

        if not model_loaded:
            self.log(
                'Pre-trained models not found — DTSS will use docking score + rules only',
                'WARN'
            )

        # DG component
        best_score = self._best_docking_score(docking_results)
        dg_norm    = self._normalise_dg(best_score)

        # ADMET component
        admet_score = self._compute_admet_score(admet)

        # Pocket druggability component
        pocket_drugg = self._best_pocket_drugg(docking_results)

        # Literature component (0–1, already normalised by chembl client)
        lit = min(max(float(lit_score), 0.0), 1.0)

        # DTSS — use reduced ML weight for data-scarce targets
        model_name = self._resolve_model_name(target_id or '')
        W = (DTSS_WEIGHTS_SPARSE
             if model_name in SPARSE_TARGETS
             else DTSS_WEIGHTS)
        # ── Standard DTSS (5 components) ─────────────────────────────
        dtss_base = (W['dg']     * dg_norm
                   + W['admet']  * admet_score
                   + W['ml']     * ml_confidence
                   + W['pocket'] * pocket_drugg
                   + W['lit']    * lit)

        # ── RNA Displacement Score (6th component, TDP-43 targets) ───
        # For RNA-binding proteins, add RNA cleft competition score.
        # Weights rescaled: base×0.90 + RNA×0.10 (weights still sum=1.0)
        rna_score = 0.5   # neutral default
        rna_mode  = False
        resolved_name = self._resolve_model_name(target_id or '')
        if (resolved_name in RNA_BINDING_TARGETS or
                str(target_id).upper() in RNA_SITE_TARGETS):
            try:
                from saidock.modules.rna_docking import RNADockingMode
                rna_docker = RNADockingMode(log_fn=self.log)
                rna_result = rna_docker.score_rna_displacement_smiles_only(
                    smiles, target_id)
                rna_score = rna_result.get('rna_displacement_score', 0.5)
                rna_mode  = rna_result.get('rna_mode', False)
            except Exception as e:
                self.log(f'RNA scoring skipped: {e}', 'WARN')

        if rna_mode:
            dtss = round(min(max(0.90 * dtss_base + 0.10 * rna_score, 0.0), 1.0), 4)
        else:
            dtss = round(min(max(dtss_base, 0.0), 1.0), 4)

        category = self._categorise(dtss)

        result = {
            'DTSS':                  dtss,
            'binding_category':      category,
            'dg_norm':               round(dg_norm,        4),
            'admet_score':           round(admet_score,    4),
            'ml_confidence':         round(ml_confidence,  4),
            'pocket_druggability':   round(pocket_drugg,   4),
            'literature_score':      round(lit,            4),
            'best_docking_score':    round(best_score,     3),
            'model_loaded':          model_loaded,
            'rna_displacement_score': round(rna_score, 4),
            'rna_mode':               rna_mode,
            'components': {
                'DG_norm (×0.35)':               round(dg_norm*0.35,       4),
                'ADMET_score (×0.20)':           round(admet_score*0.20,   4),
                'ML_confidence (×0.20)':         round(ml_confidence*0.20, 4),
                'Pocket_druggability (×0.15)':   round(pocket_drugg*0.15,  4),
                'Literature_score (×0.10)':      round(lit*0.10,           4),
                'RNA_displacement (×0.10 if TDP-43)': round(rna_score*0.10, 4),
            },
        }
        return result

    # ── Internal helpers ──────────────────────────────────────────────────────
    # PDB ID → common target name (used for model file lookup)
    PDB_TO_NAME = {
    # TDP-43 structures
    '5HDN': 'TDP43_RRM1',  '6N37': 'TDP43_RRM1',  '2N3X': 'TDP43_RRM1',
    '2N3Y': 'TDP43_RRM1',  '4BS2': 'TDP43_RRM1',  '4IUF': 'TDP43_RRM1',
    # CK2-alpha (casein kinase 2)
    '4ACC': 'CK2a',        '3NSZ': 'CK2a',         '3PE1': 'CK2a',
    '5CSH': 'CK2a',        '6RVF': 'CK2a',         '5LDU': 'CK2a',
    # CDK6
    '7XS0': 'CDK6',        '2EUF': 'CDK6',         '5L2T': 'CDK6',
    '6EZ8': 'CDK6',        '4EZ5': 'CDK6',         '1XO2': 'CDK6',
    # Keap1
    '4XMB': 'Keap1',       '4L7B': 'Keap1',        '3WN7': 'Keap1',
    '2FLU': 'Keap1',       '6CCX': 'Keap1',
    # Parkin
    '5C1Z': 'Parkin_UBL',  '4BM9': 'Parkin_UBL',   '5N2W': 'Parkin_UBL',
    # CYP2C9
    '2C9V': 'TDP43_RRM1',  # fallback — no CYP2C9 model yet
    # CZ8 series
    '6CZ8': 'TDP43_RRM1',  # fallback
}

    def _resolve_model_name(self, target_id: str):
        """Map PDB ID or model name string to _rf.pkl filename stem."""
        PDB_TO_NAME = {
            "5HDN":"TDP43_RRM1", "6N37":"TDP43_RRM1", "2N3X":"TDP43_RRM1",
            "2N3Y":"TDP43_RRM1", "4IUF":"TDP43_RRM1", "4BS2":"TDP43_RRM1",
            "2KFN":"TDP43_RRM1", "6CFH":"TDP43_RRM1", "2C9V":"TDP43_RRM1",
            "6CZ8":"TDP43_RRM1",
            "6TI9":"TDP43_pept", "5WHN":"TDP43_pept",
            "4ACC":"CK2a", "3NSZ":"CK2a", "3PE1":"CK2a", "6RVF":"CK2a",
            "3KXM":"CK2a", "4MD7":"CK2a", "6HBT":"CK2a", "3OWJ":"CK2a",
            "2PVR":"CK2a", "5LDU":"CK2a",
            "7XS0":"CDK6",  "6EZ8":"CDK6",  "5L2I":"CDK6",  "2EUF":"CDK6",
            "1XO2":"CDK6",  "2NHD":"CDK6",  "4EZ5":"CDK6",
            "4XMB":"Keap1", "4L7B":"Keap1", "5FNQ":"Keap1", "6QIB":"Keap1",
            "4ZY3":"Keap1", "5CGJ":"Keap1",
            "5C1Z":"Parkin_UBL", "4BM9":"Parkin_UBL", "4K7D":"Parkin_UBL",
            "5N2W":"Parkin_UBL",
        }
        tid = str(target_id).upper().strip()
        # 1. PDB lookup
        name = PDB_TO_NAME.get(tid)
        if name and (self.model_dir / f"{name}_rf.pkl").exists():
            return name
        # 2. Direct model name (target_id IS already the model stem)
        if (self.model_dir / f"{target_id}_rf.pkl").exists():
            return target_id
        # 3. Case-insensitive partial match
        for c in self.model_dir.glob("*_rf.pkl"):
            stem = c.stem.replace("_rf", "")
            if stem.lower() in tid.lower() or tid.lower() in stem.lower():
                return stem
        # 4. Fuzzy fallback — any available model
        cands = list(self.model_dir.glob("*_rf.pkl"))
        if cands:
            best = cands[0].stem.replace("_rf", "")
            self.log(f"No model for {target_id} — using {best}", "WARN")
            return best
        return None

    def _applicability_domain(self, smiles: str, model_name: str) -> str:
        """Tanimoto-based applicability domain check."""
        try:
            from rdkit import Chem, DataStructs
            from rdkit.Chem import AllChem
            import numpy as np
            mol = Chem.MolFromSmiles(smiles)
            if mol is None: return 'UNKNOWN'
            fp_query = AllChem.GetMorganFingerprintAsBitVect(mol, 2, 512)

            # Load training SMILES for this model if available
            train_file = self.model_dir / f'{model_name}_train_smiles.json'
            if not train_file.exists():
                return 'IN_DOMAIN'   # cannot check → assume in domain

            import json
            train_smiles = json.loads(train_file.read_text())
            fps_train = []
            for s in train_smiles[:200]:
                m = Chem.MolFromSmiles(s)
                if m: fps_train.append(
                    AllChem.GetMorganFingerprintAsBitVect(m, 2, 512))
            if not fps_train: return 'IN_DOMAIN'

            sims = DataStructs.BulkTanimotoSimilarity(fp_query, fps_train)
            max_sim = max(sims)
            # Threshold: IN_DOMAIN if max similarity to any training compound ≥ 0.30
            return 'IN_DOMAIN' if max_sim >= 0.30 else 'OUT_OF_DOMAIN'
        except Exception:
            return 'UNKNOWN'

    def _predict(self, smiles: str, target_id: str):
        """Predict P(active) using 524-feature calibrated RF model.
        Feature vector: 512 Morgan FP (r=2) + 12 RDKit descriptors = 524.
        """
        model_name = self._resolve_model_name(target_id)
        if model_name is None:
            self.log('Pre-trained model not found — using docking score only', 'WARN')
            return 0.50, False

        model_path  = self.model_dir / f'{model_name}_rf.pkl'
        scaler_path = self.model_dir / f'{model_name}_scaler.pkl'

        if not model_path.exists():
            self.log(f'Model file missing: {model_path.name}', 'WARN')
            return 0.50, False

        try:
            import joblib
            import numpy as np
            from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors

            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                self.log(f'Invalid SMILES passed to _predict', 'WARN')
                return 0.50, False

            # 512 Morgan FP + 12 descriptors = 524 (matches training exactly)
            fp   = list(AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=512))
            desc = [
                Descriptors.MolWt(mol),
                Descriptors.MolLogP(mol),
                float(rdMolDescriptors.CalcNumHBD(mol)),
                float(rdMolDescriptors.CalcNumHBA(mol)),
                float(rdMolDescriptors.CalcTPSA(mol)),
                float(rdMolDescriptors.CalcNumRotatableBonds(mol)),
                float(Descriptors.qed(mol)),
                float(rdMolDescriptors.CalcNumAromaticRings(mol)),
                float(Descriptors.FractionCSP3(mol)),
                float(mol.GetNumHeavyAtoms()),
                float(Descriptors.MolMR(mol)),
                float(rdMolDescriptors.CalcNumRings(mol)),
            ]
            feats = np.array(fp + desc, dtype=np.float64).reshape(1, -1)

            rf     = joblib.load(model_path)
            scaler = joblib.load(scaler_path) if scaler_path.exists() else None
            if scaler is not None:
                feats = scaler.transform(feats)

            prob = float(rf.predict_proba(feats)[0][1])
            self.log(f'ML [{model_name}]: conf={prob:.4f}', 'INFO')
            return round(prob, 4), True

        except Exception as e:
            self.log(f'_predict [{model_name}] failed: {e}', 'WARN')
            return 0.50, False

    def _featurise_activities(self, activities: list):
        X, y, smiles_list = [], [], []
        for act in activities:
            smi = act.get('smiles', '')
            pchembl = act.get('pchembl_value')
            if not smi or pchembl is None:
                continue
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue
            try:
                feats = self._mol_features(mol)
            except Exception:
                continue

            pchembl = float(pchembl)
            if pchembl >= ACTIVE_THRESHOLD:
                label = 1
            elif pchembl <= INACTIVE_THRESHOLD:
                label = 0
            else:
                continue   # ambiguous zone — exclude

            X.append(feats)
            y.append(label)
            smiles_list.append(smi)

        return np.array(X), np.array(y), smiles_list

    @staticmethod
    def _mol_features(mol) -> np.ndarray:
        """Morgan FP-2048 (rdFingerprintGenerator) + 15 RDKit descriptors."""
        try:
            from rdkit.Chem import rdFingerprintGenerator
            gen = rdFingerprintGenerator.GetMorganGenerator(
                radius=FP_RADIUS, fpSize=FP_BITS)
            arr = np.array(gen.GetFingerprintAsNumPy(mol), dtype=np.float32)
        except (ImportError, AttributeError):
            fp  = AllChem.GetMorganFingerprintAsBitVect(mol, FP_RADIUS, FP_BITS)
            arr = np.zeros(FP_BITS, dtype=np.float32)
            DataStructs.ConvertToNumpyArray(fp, arr)

        desc = np.array([
            Descriptors.MolWt(mol),
            Descriptors.MolLogP(mol),
            Descriptors.NumHDonors(mol),
            Descriptors.NumHAcceptors(mol),
            Descriptors.TPSA(mol),
            rdMolDescriptors.CalcNumRotatableBonds(mol),
            rdMolDescriptors.CalcNumRings(mol),
            rdMolDescriptors.CalcNumAromaticRings(mol),
            Descriptors.FractionCSP3(mol),
            Descriptors.BertzCT(mol),
            Descriptors.Chi0v(mol),
            Descriptors.Chi1v(mol),
            Descriptors.HallKierAlpha(mol),
            Descriptors.Kappa1(mol),
            Descriptors.Kappa2(mol),
        ], dtype=np.float32)

        return np.concatenate([arr, desc])

    @staticmethod
    def _descriptor_names():
        fp_names = [f'FP_{i}' for i in range(FP_BITS)]
        desc_names = [
            'MW','LogP','HBD','HBA','TPSA','RotBonds','Rings',
            'AromaticRings','FractionCSP3','BertzCT',
            'Chi0v','Chi1v','HallKierAlpha','Kappa1','Kappa2',
        ]
        return fp_names + desc_names

    @staticmethod
    def _best_docking_score(docking_results: list) -> float:
        scores = [r.get('best_score', -4.0) for r in (docking_results or [])]
        valid  = [s for s in scores if s < 0]
        return min(valid) if valid else -4.0

    @staticmethod
    def normalise_dg(dg: float) -> float:
        return round(
            min(max((dg - DG_WORST) / (DG_BEST - DG_WORST), 0.0), 1.0), 4
        )

    @classmethod
    def _normalise_dg(cls, dg: float) -> float:
        return cls.normalise_dg(dg)

    @staticmethod
    def _compute_admet_score(admet: dict) -> float:
        if not admet:
            return 0.5
        score = 0.0
        # Lipinski (25%)
        viol = admet.get('Lipinski_violations', 0) or 0
        score += 0.25 * max(0.0, 1.0 - viol * 0.5)
        # QED (20%)
        qed = float(admet.get('QED', 0.5) or 0.5)
        score += 0.20 * qed
        # GI absorption (20%)
        gi = admet.get('GI_absorption', 'Low') or 'Low'
        score += 0.20 * (1.0 if 'high' in str(gi).lower() else 0.3)
        # TPSA (15%) — ideal 40–90 Å²
        tpsa = float(admet.get('TPSA', 80) or 80)
        tpsa_s = 1.0 if 40 <= tpsa <= 90 else max(0.0, 1.0 - abs(tpsa - 65)/100)
        score += 0.15 * tpsa_s
        # BBB (10%)
        bbb = admet.get('BBB_penetrant', False)
        score += 0.10 * (1.0 if bbb else 0.3)
        # Solubility (10%)
        sol = str(admet.get('Solubility_class', '') or '').lower()
        sol_s = 1.0 if 'highly' in sol else (0.7 if 'moderate' in sol else 0.3)
        score += 0.10 * sol_s
        return round(min(score, 1.0), 4)

    @staticmethod
    def _best_pocket_drugg(docking_results: list) -> float:
        if not docking_results:
            return 0.5
        # Pocket with best (most negative) docking score
        valid = [r for r in docking_results if r.get('best_score', 0) < 0]
        if not valid:
            return 0.5
        best_r = min(valid, key=lambda r: r['best_score'])
        return float(best_r.get('pocket_drugg_score', 0.5) or 0.5)

    @staticmethod
    def _categorise(dtss: float) -> str:
        if   dtss >= 0.80: return 'Excellent (≥0.80)'
        elif dtss >= 0.65: return 'Good (0.65–0.80)'
        elif dtss >= 0.50: return 'Moderate (0.50–0.65)'
        elif dtss >= 0.35: return 'Weak (0.35–0.50)'
        else:              return 'Poor (<0.35)'

    # ── Added methods ────────────────────────────────────────────────
    def _compute_features(self, mol):
        """Compute 12 RDKit descriptors for a molecule."""
        try:
            from rdkit.Chem import Descriptors, rdMolDescriptors
            return [
                Descriptors.MolWt(mol),
                Descriptors.MolLogP(mol),
                Descriptors.NumHDonors(mol),
                Descriptors.NumHAcceptors(mol),
                Descriptors.TPSA(mol),
                Descriptors.NumRotatableBonds(mol),
                rdMolDescriptors.CalcNumRings(mol),
                rdMolDescriptors.CalcNumAromaticRings(mol),
                Descriptors.FractionCSP3(mol),
                Descriptors.NumAliphaticRings(mol),
                Descriptors.HeavyAtomCount(mol),
                Descriptors.NumHeteroatoms(mol),
            ]
        except Exception:
            return None

    def _scaffold_split_eval(self, X, y, valid_smiles: list) -> float:
        """Bemis-Murcko scaffold split AUC. Falls back to random split."""
        import numpy as np
        from sklearn.ensemble import RandomForestClassifier
        from sklearn.metrics import roc_auc_score
        try:
            from rdkit import Chem
            from rdkit.Chem.Scaffolds import MurckoScaffold
            scaffolds = {}
            for i, smi in enumerate(valid_smiles):
                mol = Chem.MolFromSmiles(smi)
                scaf = (MurckoScaffold.MurckoScaffoldSmiles(
                            mol=mol, includeChirality=False)
                        if mol else "")
                scaffolds.setdefault(scaf, []).append(i)
            sets = sorted(scaffolds.values(), key=len, reverse=True)
            test_idx, train_idx = [], []
            tgt = max(2, int(len(y) * 0.2))
            for ss in sets:
                (test_idx if len(test_idx) < tgt else train_idx).extend(ss)
            if len(test_idx) < 2 or len(np.unique(y[test_idx])) < 2:
                raise ValueError("not enough scaffold diversity")
            rf = RandomForestClassifier(100, class_weight="balanced",
                                        n_jobs=-1, random_state=42)
            rf.fit(X[train_idx], y[train_idx])
            auc = roc_auc_score(y[test_idx],
                                rf.predict_proba(X[test_idx])[:,1])
            self.log(f"  Scaffold AUC={auc:.3f} "
                     f"train={len(train_idx)} test={len(test_idx)}")
            return float(auc)
        except Exception as e:
            self.log(f"  Scaffold eval fallback: {e}", "WARN")
            try:
                from sklearn.model_selection import train_test_split
                Xtr,Xte,ytr,yte = train_test_split(
                    X, y, test_size=0.25, stratify=y, random_state=42)
                rf2 = RandomForestClassifier(100, class_weight="balanced",
                                             n_jobs=-1, random_state=42)
                rf2.fit(Xtr, ytr)
                return float(roc_auc_score(yte, rf2.predict_proba(Xte)[:,1]))
            except Exception:
                return 0.5

    def predict(self, smiles: str, target_id: str) -> float:
        """Public alias returning activity probability (0–1)."""
        prob, _ = self._predict(smiles, target_id)
        return prob