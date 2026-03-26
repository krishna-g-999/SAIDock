#!/usr/bin/env python3
"""
SAIDock Batch Pipeline
Modes:
  screen : 1 target × N ligands (virtual screening)
  scan   : 1 ligand × N targets (polypharmacology/off-target)
  matrix : N ligands × M targets (cross-screening)
"""
import os, sys, json, time, csv, traceback
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field, asdict
from typing import List, Dict, Optional
import argparse


# ── Result container ─────────────────────────────────────────────────────────
@dataclass
class DockingResult:
    compound_name: str
    smiles:        str
    target:        str
    best_score:    float  = 0.0
    dtss:          float  = 0.0
    category:      str    = ''
    ml_confidence: float  = 0.0
    admet_qed:     float  = 0.0
    admet_bbb:     bool   = False
    admet_gi:      str    = ''
    hbonds:        int    = 0
    n_poses:       int    = 0
    report_html:   str    = ''
    error:         str    = ''
    elapsed_s:     float  = 0.0


# ── Single run wrapper (runs in subprocess) ───────────────────────────────────
def _run_single(task: dict) -> dict:
    """Callable for ProcessPoolExecutor — runs one saidock pipeline."""
    import argparse, time
    t0 = time.time()
    try:
        from saidock.pipeline import SAIDockRun

        # Build a minimal args namespace from the task dict
        args = argparse.Namespace(
            smiles          = task['smiles'],
            name            = task['name'],
            target          = task['target'],
            output          = task['output'],
            chain           = task.get('chain', 'A'),
            exhaustiveness  = task.get('exhaustiveness', 16),
            n_pockets       = task.get('n_pockets', 3),
            cpu             = task.get('cpu', 2),
            engine          = task.get('engine', 'vina'),
            no_interactions = task.get('no_interactions', False),
            no_admet        = False,
            no_ml           = False,
        )

        run = SAIDockRun(args)
        run.execute()
        state = run.state

        dr = DockingResult(
            compound_name = task['name'],
            smiles        = task['smiles'],
            target        = task['target'],
            best_score    = float(state.get('ml_results', {}).get('best_docking_score', 0)),
            dtss          = float(state.get('ml_results', {}).get('DTSS', 0)),
            category      = state.get('ml_results', {}).get('binding_category', ''),
            ml_confidence = float(state.get('ml_results', {}).get('ml_confidence', 0)),
            admet_qed     = float((state.get('admet') or {}).get('QED', 0)),
            admet_bbb     = bool((state.get('admet') or {}).get('BBB_penetrant', False)),
            admet_gi      = str((state.get('admet') or {}).get('GI_absorption', '')),
            hbonds        = int((state.get('interactions') or {}).get('hbonds', 0)),
            n_poses       = sum(r.get('n_poses', 0)
                                for r in (state.get('docking_results') or [])),
            report_html   = str(Path(task['output']) / 'report' / 'saidock_report.html'),
            elapsed_s     = round(time.time() - t0, 1),
        )
        return asdict(dr)

    except Exception as e:
        return asdict(DockingResult(
            compound_name=task.get('name','?'),
            smiles=task.get('smiles',''),
            target=task.get('target','?'),
            error=str(e),
            elapsed_s=round(time.time()-t0, 1),
        ))


# ── Main Batch orchestrator ───────────────────────────────────────────────────
class BatchPipeline:
    def __init__(self, args):
        self.args    = args
        self.outdir  = Path(args.output).resolve()
        self.outdir.mkdir(parents=True, exist_ok=True)
        self.results: List[DockingResult] = []

    def run(self):
        mode = self.args.mode
        print(f'\nSAIDock Batch | mode={mode} | workers={self.args.workers}')
        print('─' * 60)

        if   mode == 'screen': tasks = self._build_screen_tasks()
        elif mode == 'scan':   tasks = self._build_scan_tasks()
        elif mode == 'matrix': tasks = self._build_matrix_tasks()
        else:
            print(f'Unknown mode: {mode}'); return

        print(f'Tasks queued: {len(tasks)}')
        self._execute(tasks)
        self._save_results()
        self._print_summary()

    # ── Task builders ─────────────────────────────────────────────────────────
    def _build_screen_tasks(self) -> list:
        """1 target × N ligands."""
        if not self.args.library:
            raise ValueError('--library required for screen mode')
        if not self.args.target:
            raise ValueError('--target required for screen mode')
        compounds = self._load_library(self.args.library)
        target    = self.args.target
        return [self._task(c['smiles'], c['name'], target,
                           f"{c['name'].replace(' ','_')}_{target}")
                for c in compounds]

    def _build_scan_tasks(self) -> list:
        """1 ligand × N targets."""
        if not self.args.smiles:
            raise ValueError('--smiles required for scan mode')
        targets = self._load_targets(self.args.targets)
        name    = self.args.name
        return [self._task(self.args.smiles, name, t,
                           f"{name.replace(' ','_')}_{t}")
                for t in targets]

    def _build_matrix_tasks(self) -> list:
        """N ligands × M targets."""
        compounds = self._load_library(self.args.library)
        targets   = self._load_targets(self.args.targets)
        tasks = []
        for c in compounds:
            for t in targets:
                tasks.append(self._task(
                    c['smiles'], c['name'], t,
                    f"{c['name'].replace(' ','_')}_{t}"
                ))
        return tasks

    def _task(self, smiles, name, target, run_id) -> dict:
        return {
            'smiles':         smiles,
            'name':           name,
            'target':         target,
            'output':         str(self.outdir / run_id),
            'exhaustiveness': self.args.exhaustiveness,
            'n_pockets':      self.args.n_pockets,
            'cpu':            max(1, self.args.cpu // self.args.workers),
            'engine':         self.args.engine,
        }

    # ── Execution ─────────────────────────────────────────────────────────────
    def _execute(self, tasks: list):
        t0      = time.time()
        n       = len(tasks)
        done    = 0
        workers = min(self.args.workers, n)

        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = {pool.submit(_run_single, t): t for t in tasks}
            for future in as_completed(futures):
                done += 1
                result = future.result()
                self.results.append(result)
                elapsed = time.time() - t0
                eta     = (elapsed / done) * (n - done)
                status  = 'ERROR' if result.get('error') else 'OK'
                print(f"  [{done:>3}/{n}] [{status}] "
                      f"{result['compound_name']:<18} vs {result['target']:<6} "
                      f"DTSS={result['dtss']:.3f}  "
                      f"Vina={result['best_score']:.2f}  "
                      f"ETA {eta/60:.1f}min")

    # ── Output ────────────────────────────────────────────────────────────────
    def _save_results(self):
        # Sort by DTSS descending
        self.results.sort(key=lambda r: r.get('dtss', 0), reverse=True)

        # JSON
        json_out = self.outdir / 'batch_results.json'
        with open(json_out, 'w') as f:
            json.dump(self.results, f, indent=2)

        # CSV
        csv_out = self.outdir / 'batch_results.csv'
        if self.results:
            keys = list(self.results[0].keys())
            with open(csv_out, 'w', newline='') as f:
                w = csv.DictWriter(f, fieldnames=keys)
                w.writeheader()
                w.writerows(self.results)

        print(f'\n  Results saved:')
        print(f'    {csv_out}')
        print(f'    {json_out}')

        # HTML summary table
        self._write_summary_html()

    def _write_summary_html(self):
        html_out = self.outdir / 'batch_summary.html'
        rows = ''
        for r in self.results:
            color = ('#d4edda' if r.get('dtss', 0) >= 0.65
                     else '#fff3cd' if r.get('dtss', 0) >= 0.50
                     else '#f8d7da')
            err   = f'<br><small style="color:red">{r["error"]}</small>' if r.get('error') else ''
            rep   = (f'<a href="{r["report_html"]}">📄</a>'
                     if r.get('report_html') and not r.get('error') else '—')
            rows += f"""<tr style="background:{color}">
  <td>{r['compound_name']}</td><td>{r['target']}</td>
  <td>{r['best_score']:.3f}</td><td><b>{r['dtss']:.3f}</b></td>
  <td>{r.get('category','')}</td><td>{r.get('ml_confidence',0):.3f}</td>
  <td>{r.get('admet_qed',0):.3f}</td>
  <td>{'✅' if r.get('admet_bbb') else '❌'}</td>
  <td>{r.get('admet_gi','')}</td>
  <td>{r.get('hbonds',0)}</td>
  <td>{r.get('elapsed_s',0):.0f}s</td>
  <td>{rep}{err}</td>
</tr>"""

        n_good = sum(1 for r in self.results if r.get('dtss', 0) >= 0.65)
        html = f"""<!DOCTYPE html><html><head>
<title>SAIDock Batch Results</title>
<style>
  body{{font-family:monospace;padding:20px}}
  table{{border-collapse:collapse;width:100%;font-size:13px}}
  th{{background:#343a40;color:white;padding:8px;position:sticky;top:0}}
  td{{padding:6px;border-bottom:1px solid #dee2e6}}
  h2{{color:#343a40}}
</style></head><body>
<h2>SAIDock Batch Results — {self.args.mode.upper()} mode</h2>
<p><b>{len(self.results)}</b> runs completed &nbsp;|&nbsp;
   <b>{n_good}</b> hits (DTSS ≥ 0.65)</p>
<table><thead><tr>
  <th>Compound</th><th>Target</th><th>Vina (kcal/mol)</th><th>DTSS</th>
  <th>Category</th><th>ML Conf.</th><th>QED</th><th>BBB</th><th>GI</th>
  <th>H-bonds</th><th>Time</th><th>Report</th>
</tr></thead><tbody>{rows}</tbody></table>
</body></html>"""
        html_out.write_text(html)
        print(f'    {html_out}')

    def _print_summary(self):
        n   = len(self.results)
        ok  = [r for r in self.results if not r.get('error')]
        top = [r for r in ok if r.get('dtss', 0) >= 0.65]
        print(f'\n{"="*60}')
        print(f'BATCH COMPLETE: {len(ok)}/{n} OK  |  {len(top)} hits (DTSS≥0.65)')
        if top:
            print('\nTop 5 hits:')
            for r in top[:5]:
                print(f"  {r['compound_name']:<20} vs {r['target']:<6} "
                      f"DTSS={r['dtss']:.3f}  Vina={r['best_score']:.3f}  "
                      f"ML={r['ml_confidence']:.3f}")
        print(f'{"="*60}')

    # ── Helpers ───────────────────────────────────────────────────────────────
    @staticmethod
    def _load_library(path: str) -> List[Dict]:
        """Load CSV with columns: smiles, name."""
        compounds = []
        with open(path, newline='') as f:
            # Auto-detect delimiter
            sample  = f.read(2048); f.seek(0)
            dialect = csv.Sniffer().sniff(sample, delimiters=',\t;')
            reader  = csv.DictReader(f, dialect=dialect)
            for row in reader:
                row = {k.strip().lower(): v.strip() for k, v in row.items()}
                smiles = row.get('smiles') or row.get('smi') or ''
                name   = (row.get('name') or row.get('compound_name')
                          or row.get('id') or 'compound')
                if smiles:
                    compounds.append({'smiles': smiles, 'name': name})
        return compounds

    @staticmethod
    def _load_targets(targets_arg) -> List[str]:
        """Accept list of PDB IDs or path to a text file."""
        if not targets_arg:
            raise ValueError('No targets provided')
        result = []
        for t in targets_arg:
            p = Path(t)
            if p.exists() and p.suffix in ('.txt', '.csv'):
                result.extend(line.strip() for line in p.read_text().splitlines()
                               if line.strip() and not line.startswith('#'))
            else:
                result.append(t.upper())
        return result
