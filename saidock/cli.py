#!/usr/bin/env python3
"""SAIDock command-line interface — all subcommands."""
import argparse, sys

def build_parser():
    parser = argparse.ArgumentParser(
        prog='saidock',
        description='SAIDock v1.0.0 — Automated Drug-Target Docking & Assessment',
    )
    sub = parser.add_subparsers(dest='command')

    # ── saidock run ──────────────────────────────────────────────────────────
    run = sub.add_parser('run', help='Single compound vs single target')
    run.add_argument('--smiles',          required=True)
    run.add_argument('--name',            default='compound')
    run.add_argument('--target',          required=True)
    run.add_argument('--output',          required=True)
    run.add_argument('--chain',           default='A')
    run.add_argument('--exhaustiveness',  type=int,   default=32)
    run.add_argument('--n-pockets',       type=int,   default=5)
    run.add_argument('--cpu',             type=int,   default=4)
    run.add_argument('--engine',          default='vina',
                     choices=['vina','gnina','both'],
                     help='Docking engine: vina (fast), gnina (accurate), both (tiered)')
    run.add_argument('--no-interactions', action='store_true')
    run.add_argument('--no-admet',        action='store_true')
    run.add_argument('--no-ml',           action='store_true')

    # ── saidock batch ─────────────────────────────────────────────────────────
    batch = sub.add_parser('batch', help='Batch docking (screen / scan / matrix)')
    batch.add_argument('--mode', required=True,
                       choices=['screen','scan','matrix'],
                       help='screen=1 target×N ligands | scan=1 ligand×N targets | matrix=N×M')
    # screen / matrix
    batch.add_argument('--library',   default=None,
                       help='CSV with columns: smiles,name  (for screen/matrix modes)')
    # scan / matrix
    batch.add_argument('--targets',   default=None, nargs='+',
                       help='PDB IDs or path to targets.txt  (for scan/matrix modes)')
    # scan: single ligand
    batch.add_argument('--smiles',    default=None)
    batch.add_argument('--name',      default='compound')
    # shared
    batch.add_argument('--target',    default=None,
                       help='Single target PDB ID (for screen mode)')
    batch.add_argument('--output',    required=True)
    batch.add_argument('--exhaustiveness', type=int, default=16)
    batch.add_argument('--n-pockets', type=int, default=3)
    batch.add_argument('--cpu',       type=int, default=4)
    batch.add_argument('--workers',   type=int, default=2,
                       help='Parallel pipeline workers')
    batch.add_argument('--engine',    default='vina',
                       choices=['vina','gnina','both'])
    batch.add_argument('--top-n',     type=int, default=0,
                       help='Re-dock top N hits with gnina (when --engine both)')

    # ── saidock surface ───────────────────────────────────────────────────────
    surface = sub.add_parser('surface',
                             help='Pan-pocket surface analysis & functional annotation')
    surface.add_argument('--target',  required=True,
                         help='PDB ID or path to local PDB file')
    surface.add_argument('--output',  required=True)
    surface.add_argument('--chain',   default='A')
    surface.add_argument('--uniprot', default=None,
                         help='UniProt accession for functional annotation (auto-detected if omitted)')
    surface.add_argument('--min-vol', type=float, default=100.0,
                         help='Minimum pocket volume in Å³ to include (default 100)')
    surface.add_argument('--max-pockets', type=int, default=20)

    # ── saidock train ─────────────────────────────────────────────────────────
    train = sub.add_parser('train',
                           help='(Re-)train ML models for a target from ChEMBL data')
    train.add_argument('--target',        required=True,
                       help='ChEMBL target ID (e.g. CHEMBL2185) or PDB ID')
    train.add_argument('--name',          required=True,
                       help='Model name used for saving (e.g. CK2a)')
    train.add_argument('--model-dir',     default='models/')
    train.add_argument('--data-file', default=None,
                   help='CSV with smiles,pchembl_value columns ')
    train.add_argument('--max-records',   type=int, default=1000)
    train.add_argument('--eval',          action='store_true',
                       help='Report AUPRC + MCC + scaffold-split AUC after training')

    return parser


def main():
    parser = build_parser()
    args   = parser.parse_args()

    if not args.command:
        parser.print_help()
        sys.exit(0)

    if args.command == 'run':
        from saidock.pipeline import SAIDockRun
        SAIDockRun(args).execute()

    elif args.command == 'batch':
        from saidock.batch_pipeline import BatchPipeline
        BatchPipeline(args).run()

    elif args.command == 'surface':
        from saidock.surface_pipeline import SurfacePipeline
        SurfacePipeline(args).run()

    elif args.command == 'train':
        from saidock.modules.ml_scoring import MLScorer
        from saidock.utils.chembl_client import fetch_target_activities
        print(f'Fetching ChEMBL data for {args.target}...')
        # --target is optional when --data-file is provided
        if args.data_file and not args.target:
            args.target = args.name   # use model name as target label
        data = fetch_target_activities(args.target, max_records=args.max_records) if args.target else []
        if getattr(args, 'data_file', None):
            import csv as _csv
            with open(args.data_file) as _f:
                extra = [{'smiles': r['smiles'],
                           'pchembl_value': float(r['pchembl_value'])}
                          for r in _csv.DictReader(_f)]
            data = (data or []) + extra
            print(f'  + {len(extra)} records from {args.data_file}')
        print(f'  {len(data)} activities fetched')
        scorer  = MLScorer(model_dir=args.model_dir)
        metrics = scorer.train(data, target_id=args.name)
        import json
        print(json.dumps(metrics, indent=2))


if __name__ == '__main__':
    main()
