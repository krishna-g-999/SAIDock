"""Entry point for `saidock surface` command."""
from pathlib import Path
from saidock.utils.logger import SAILogger
from saidock.modules.surface_analyser import SurfaceAnalyser


class SurfacePipeline:
    def __init__(self, args):
        self.args   = args
        self.outdir = Path(args.output).resolve()
        self.outdir.mkdir(parents=True, exist_ok=True)
        self.logger = SAILogger(
            log_path=str(self.outdir / 'surface_analysis.log'),
            verbose=True,
        )

    def run(self):
        self.logger.step(f'SAIDock Surface Analysis | target={self.args.target}')

        # 1. Prepare receptor PDB
        self.logger.step('[1/3] Target Preparation')
        from saidock.modules.target_prep import TargetPrep
        tp = TargetPrep(
            target = self.args.target,
            chain  = getattr(self.args, 'chain', 'A'),
            outdir = str(self.outdir / 'target'),
            logger = self.logger,
        )
        receptor_pdb = str(tp.prepare())

        # 2. Run surface analysis
        self.logger.step('[2/3] Pan-Pocket Surface Analysis')
        sa = SurfaceAnalyser(
            receptor_pdb  = receptor_pdb,
            target_id     = tp.target_id,
            uniprot_id    = getattr(self.args, 'uniprot', None),
            outdir        = str(self.outdir),
            min_vol       = getattr(self.args, 'min_vol', 100.0),
            max_pockets   = getattr(self.args, 'max_pockets', 20),
            logger        = self.logger,
        )
        atlas = sa.analyse()

        # 3. Report
        self.logger.step('[3/3] Report Generation')
        sa.write_report(atlas)
        self.logger.ok(
            f'Surface analysis complete | {len(atlas)} pockets analysed | '
            f'Report: {self.outdir}/surface_report.html'
        )
