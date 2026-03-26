import os
import sys
import json
import time
from pathlib import Path
from saidock.utils.logger import SAILogger

PIPELINE_STAGES = [
    'target_preparation',
    'pocket_detection',
    'ligand_preparation',
    'docking',
    'interactions',
    'admet',
    'ml_dtss',
    'report',
]

class SAIDockRun:
    """
    Master orchestrator for a single SAIDock assessment.
    Calls modules 01-08 in sequence. On failure of non-critical
    stages, continues with graceful degradation.
    """

    CRITICAL_STAGES = {'target_preparation', 'ligand_preparation', 'docking'}

    def __init__(self, args):
        self.args    = args
        self.outdir  = Path(args.output).resolve()
        self.outdir.mkdir(parents=True, exist_ok=True)
        for subdir in ['target','ligand','pockets','docking',
                       'interactions','admet','ml','figures','report','logs']:
            (self.outdir / subdir).mkdir(exist_ok=True)
        self.logger  = SAILogger(
            log_path = str(self.outdir / 'logs' / 'saidock_run.log'),
            verbose  = True,
        )
        self.state   = {'outdir': str(self.outdir)}
        self.t_start = time.time()

    def execute(self):
        self.logger.step(f'SAIDock v1.0.0 | Output: {self.outdir}')
        stages = [
            ('Target Preparation',  self._stage_target_prep),
            ('Pocket Detection',    self._stage_pocket_detect),
            ('Ligand Preparation',  self._stage_ligand_prep),
            ('Docking',             self._stage_docking),
            ('Interactions',        self._stage_interactions),
            ('ADMET',               self._stage_admet),
            ('ML + DTSS',           self._stage_ml_dtss),
            ('Report Generation',   self._stage_report),
        ]
        for i, (name, func) in enumerate(stages, 1):
            self.logger.step(f'[{i}/{len(stages)}] {name}')
            try:
                func()
            except Exception as exc:
                self.logger.error(f'{name} failed: {exc}')
                import traceback
                with open(self.outdir / 'logs' / 'saidock_run.log', 'a') as f:
                    traceback.print_exc(file=f)
                key = name.lower().replace(' ', '_').replace('+_','')
                if key in self.CRITICAL_STAGES:
                    self.logger.error('Critical stage failed. Aborting.')
                    sys.exit(1)
                self.logger.warn(f'Non-critical stage skipped: {name}')

        elapsed = time.time() - self.t_start
        self.logger.ok(
            f'Assessment complete in {elapsed:.1f}s | '
            f'Report: {self.outdir}/report/saidock_report.html'
        )
        with open(self.outdir / 'saidock_state.json', 'w') as f:
            json.dump(
                {k: (str(v) if not isinstance(v, (int,float,str,list,dict,bool,type(None)))
                     else v)
                 for k, v in self.state.items()},
                f, indent=2
            )

    def _stage_target_prep(self):
        from saidock.modules.target_prep import TargetPrep
        tp = TargetPrep(
            target  = self.args.target,
            chain   = getattr(self.args, 'chain', 'A'),
            outdir  = str(self.outdir / 'target'),
            logger  = self.logger,
        )
        self.state['receptor_pdb']   = str(tp.prepare())
        self.state['receptor_pdbqt'] = str(tp.convert_to_pdbqt())
        self.state['target_id']      = tp.target_id
        self.state['n_residues']     = tp.n_residues

    def _stage_pocket_detect(self):
        from saidock.modules.pocket_detect import PocketDetector
        pd = PocketDetector(
            receptor_pdb = self.state['receptor_pdb'],
            outdir       = str(self.outdir / 'pockets'),
            n_pockets    = getattr(self.args, 'n_pockets', 5),
            logger       = self.logger,
        )
        self.state['pockets'] = pd.detect()

    def _stage_ligand_prep(self):
        from saidock.modules.ligand_prep import LigandPrep
        kwargs = {'outdir': str(self.outdir / 'ligand'), 'logger': self.logger}
        if getattr(self.args, 'ligand', None):  kwargs['name']     = self.args.ligand
        if getattr(self.args, 'smiles', None):  kwargs['smiles']   = self.args.smiles
        if getattr(self.args, 'name',   None):  kwargs['name']     = self.args.name
        if getattr(self.args, 'sdf',    None):  kwargs['sdf_path'] = self.args.sdf
        lp = LigandPrep(**kwargs)
        self.state['ligand_pdbqt'] = str(lp.prepare())
        self.state['ligand_smiles']= lp.smiles
        self.state['ligand_name']  = lp.name
        self.state['ligand_mw']    = lp.mw
        # Save SMILES for PolarDock fallback
        (self.outdir / 'ligand' / 'ligand_smiles.txt').write_text(lp.smiles)

    def _stage_docking(self):
        from saidock.modules.docking import MultiPocketDocker
        docker = MultiPocketDocker(
            receptor_pdbqt  = self.state['receptor_pdbqt'],
            ligand_pdbqt    = self.state['ligand_pdbqt'],
            pockets         = self.state['pockets'],
            outdir          = str(self.outdir / 'docking'),
            exhaustiveness  = getattr(self.args, 'exhaustiveness', 32),
            cpu             = getattr(self.args, 'cpu', 4),
            logger          = self.logger,
        )
        self.state['docking_results'] = docker.run()

    def _stage_interactions(self):
        if getattr(self.args, 'no_interactions', False):
            self.logger.warn('Interaction analysis skipped (--no-interactions)')
            return
        from saidock.modules.interactions import InteractionAnalyzer
        ia = InteractionAnalyzer(
            receptor_pdb    = self.state['receptor_pdb'],
            docking_results = self.state['docking_results'],
            outdir          = str(self.outdir / 'interactions'),
            logger          = self.logger,
        )
        self.state['interactions'] = ia.analyse()

    def _stage_admet(self):
        if getattr(self.args, 'no_admet', False):
            self.logger.warn('ADMET skipped (--no-admet)')
            return
        from saidock.modules.admet import ADMETCalculator
        ac = ADMETCalculator(
            smiles = self.state['ligand_smiles'],
            name   = self.state['ligand_name'],
            outdir = str(self.outdir / 'admet'),
            logger = self.logger,
        )
        self.state['admet'] = ac.calculate()

    def _stage_ml_dtss(self):
        if getattr(self.args, 'no_ml', False):
            self.logger.warn('ML scoring skipped (--no-ml)')
            return
        from saidock.modules.ml_scoring import MLScorer
        # Fetch ChEMBL literature score
        lit_score   = 0.0
        target_str  = self.state.get('target_id', '')
        try:
            from saidock.utils.chembl_client import fetch_target_activities
            acts = fetch_target_activities(target_str, max_records=200)
            if acts:
                best_pc = max(a['pchembl_value'] for a in acts)
                self.logger.info(f'ChEMBL: best pChEMBL = {best_pc:.2f}')
                lit_score = round(min((best_pc - 4.0) / 6.0, 1.0), 4)
            else:
                self.logger.info('Not found in ChEMBL — literature_score = 0.0')
        except Exception as e:
            self.logger.warn(f'ChEMBL query failed: {e}')

        ml = MLScorer(
            ligand_smiles   = self.state['ligand_smiles'],
            target_id       = self.state['target_id'],
            docking_results = self.state['docking_results'],
            admet           = self.state.get('admet', {}),
            pockets         = self.state.get('pockets', []),
            lit_score       = lit_score,
            outdir          = str(self.outdir / 'ml'),
            logger          = self.logger,
        )
        result = ml.score()

        # Save to ml/ directory
        import json
        ml_out = self.outdir / 'ml' / 'ml_dtss_results.json'
        with open(ml_out, 'w') as f:
            json.dump(result, f, indent=2)

        self.logger.ok(
            f'DTSS = {result["DTSS"]:.3f}  Category: {result["binding_category"]}  '
            f'[DG={result["dg_norm"]:.3f} ADMET={result["admet_score"]:.3f} '
            f'ML={result["ml_confidence"]:.3f} Pocket={result["pocket_druggability"]:.3f} '
            f'Lit={result["literature_score"]:.3f}]'
        )
        self.state['ml_results'] = result

    def _stage_report(self):
        from saidock.modules.report import ReportGenerator
        rg = ReportGenerator(
            state  = self.state,
            outdir = str(self.outdir / 'report'),
            logger = self.logger,
        )
        rg.generate()
