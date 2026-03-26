"""SAIDock Web API — fixed job lifecycle, correct output paths."""
import uuid, subprocess, json, asyncio, time, threading, os
from pathlib import Path
from fastapi import FastAPI, Form, HTTPException
from fastapi.staticfiles import StaticFiles
from fastapi.responses import FileResponse, HTMLResponse, StreamingResponse
from fastapi.middleware.cors import CORSMiddleware

SAIDOCK_DIR = Path(__file__).parent.parent.parent   # ~/SAIDock
JOBS_DIR    = Path(__file__).parent.parent / 'jobs'   # web/jobs/
JOBS_DIR.mkdir(parents=True, exist_ok=True)

app = FastAPI(title="SAIDock", version="1.0.0",
              docs_url="/api/docs", redoc_url=None)
app.add_middleware(CORSMiddleware, allow_origins=["*"],
                   allow_methods=["*"], allow_headers=["*"])


# ── Background thread: run saidock and update status when done ────────
def _worker(job_id: str, smiles: str, target: str):
    job_dir = JOBS_DIR / job_id
    log     = open(job_dir / "run.log", "w", buffering=1)
    try:
        proc = subprocess.Popen(
            ["conda", "run", "--no-capture-output", "-n", "saidock",
             "saidock", "run",
             "--smiles",  smiles,
             "--target",  target,
             "--output",  str(job_dir / "results")],
            stdout=log, stderr=log,
            cwd=str(SAIDOCK_DIR)
        )
        (job_dir / "pid.txt").write_text(str(proc.pid))
        ret = proc.wait()                           # ← blocks until done

        # Detect success: report file must exist
        report = _find_report(job_dir)
        if ret == 0 and report:
            (job_dir / "status.txt").write_text("done")
        else:
            (job_dir / "status.txt").write_text("error")
    except Exception as e:
        (job_dir / "status.txt").write_text("error")
        log.write(f"\n[FATAL] Worker exception: {e}\n")
    finally:
        log.close()


def _find_report(job_dir: Path):
    """Locate saidock_report.html anywhere under job_dir."""
    if not job_dir.exists():
        return None
    # Try known locations first (fast path)
    for candidate in [
        job_dir / "results" / "report" / "saidock_report.html",
        job_dir / "results" / "saidock_report.html",
        job_dir / "saidock_report.html",
    ]:
        if candidate.exists():
            return str(candidate)
    # Fallback: recursive search
    hits = list(job_dir.rglob("saidock_report.html"))
    return str(hits[0]) if hits else None


def _find_csv(job_dir: Path) -> Path | None:
    """Find the results CSV wherever saidock wrote it."""
    candidates = [
        job_dir / "results" / "results.csv",
        job_dir / "results" / "full_results_final.csv",
        job_dir / "results" / "docking" / "results.csv",
        job_dir / "results" / "report"  / "results.csv",
    ]
    for p in job_dir.rglob("*.csv"):
        candidates.append(p)
    return next((p for p in candidates if p.exists()), None)


# ── Submit job ────────────────────────────────────────────────────────
@app.post("/api/run")
async def submit_run(smiles: str = Form(...), target: str = Form(...)):
    job_id  = uuid.uuid4().hex[:8]
    job_dir = JOBS_DIR / job_id
    job_dir.mkdir()

    target = target.upper().strip()
    params = {"smiles": smiles, "target": target, "job_id": job_id}
    (job_dir / "params.json").write_text(json.dumps(params, indent=2))
    (job_dir / "status.txt").write_text("running")

    # Fire background thread — process.wait() inside will update status
    t = threading.Thread(target=_worker, args=(job_id, smiles, target),
                         daemon=True)
    t.start()
    return {"job_id": job_id, "status": "running",
            "poll": f"/api/jobs/{job_id}"}


# ── Job status ────────────────────────────────────────────────────────
@app.get("/api/jobs/{job_id}")
def job_status(job_id: str):
    job_dir = JOBS_DIR / job_id
    if not job_dir.exists():
        raise HTTPException(404, "Job not found")

    status = (job_dir / "status.txt").read_text().strip()
    params = json.loads((job_dir / "params.json").read_text())

    # Extra safety: if log says "Assessment complete" but status wasn't set
    if status == "running":
        log_f = job_dir / "run.log"
        if log_f.exists():
            tail = log_f.read_text(errors="replace")[-500:]
            if "Assessment complete" in tail:
                (job_dir / "status.txt").write_text("done")
                status = "done"
            elif "failed. (See above for error)" in tail or \
                 "required: --output" in tail:
                (job_dir / "status.txt").write_text("error")
                status = "error"

    response = {"job_id": job_id, "status": status, **params}

    # Timing
    log_f = job_dir / "run.log"
    if log_f.exists():
        elapsed = round(time.time() - log_f.stat().st_ctime)
        response["elapsed_s"] = elapsed

    if status == "done":
        report = _find_report(job_dir)
        csv    = _find_csv(job_dir)
        response["report_url"] = f"/api/jobs/{job_id}/report" if report else None
        response["csv_url"]    = f"/api/jobs/{job_id}/csv"    if csv    else None
        if csv:
            import csv as _csv
            try:
                rows = list(_csv.DictReader(open(csv)))
                response["n_results"]  = len(rows)
                response["top_pocket"] = sorted(
                    rows, key=lambda x: float(x.get('DTSS', 0)), reverse=True
                )[0] if rows else None
            except Exception:
                response["n_results"] = 0

    return response


# ── Live log stream (SSE) ─────────────────────────────────────────────
@app.get("/api/jobs/{job_id}/log")
async def stream_log(job_id: str):
    log_path = JOBS_DIR / job_id / "run.log"
    status_f = JOBS_DIR / job_id / "status.txt"

    async def _gen():
        pos = 0
        for _ in range(1200):   # max 20 min
            if log_path.exists():
                txt = log_path.read_text(errors='replace')
                if len(txt) > pos:
                    for line in txt[pos:].splitlines():
                        yield f"data: {line}\n\n"
                    pos = len(txt)

            status = status_f.read_text().strip() if status_f.exists() else "running"

            # Catch conda CLI failure in log even before status updates
            if log_path.exists():
                tail = log_path.read_text(errors='replace')[-500:]
                if "Assessment complete" in tail:
                    status = "done"
                elif "failed. (See above for error)" in tail or \
                     "required: --output" in tail:
                    status = "error"

            if status in ("done", "error"):
                yield f"data: __JOBSTATUS__{status}\n\n"
                break
            await asyncio.sleep(1)

    return StreamingResponse(_gen(), media_type="text/event-stream",
                             headers={"Cache-Control": "no-cache",
                                      "X-Accel-Buffering": "no"})



# ── Aggregated JSON results (reads actual saidock output) ────────────
@app.get("/api/jobs/{job_id}/results")
def job_results(job_id: str):
    job_dir = JOBS_DIR / job_id
    # Primary: read saidock_state.json (complete unified result)
    state_f = job_dir / "results" / "saidock_state.json"
    if state_f.exists():
        try:
            return json.loads(state_f.read_text())
        except Exception as e:
            raise HTTPException(500, f"State parse error: {e}")

    # Fallback: assemble from individual JSON files
    rdir = job_dir / "results"
    if not rdir.exists():
        raise HTTPException(404, "Results not ready")

    result = {}
    for fname, key in [
        ("ml/ml_dtss_results.json",   "ml_results"),
        ("admet/admet_results.json",  "admet"),
        ("docking/docking_results.json", "docking_results"),
        ("pockets/pockets.json",      "pockets"),
        ("interactions/interactions.json", "interactions"),
    ]:
        fp = rdir / fname
        if fp.exists():
            try:
                result[key] = json.loads(fp.read_text())
            except:
                pass
    if not result:
        raise HTTPException(404, "No results found")
    return result


@app.get("/api/jobs/{job_id}/csv")
def job_csv(job_id: str):
    """CSV export — results embedded in full report"""
    raise HTTPException(404, "CSV export removed — use Full Report")
    return FileResponse(str(csv), media_type="text/csv",
                        filename=f"saidock_{job_id}.csv")


# ── Surface / docking HTML report ────────────────────────────────────


# ── Inline HTML report builder (no Plotly dependency) ─────────────────
@app.get("/api/jobs/{job_id}/report/{target}")
def job_report_legacy(job_id: str, target: str):
    return job_report(job_id)


# ── List all jobs ─────────────────────────────────────────────────────
@app.get("/api/jobs")
def list_jobs():
    jobs = []
    for d in sorted(JOBS_DIR.iterdir(), reverse=True)[:50]:
        sf = d / "status.txt"
        pf = d / "params.json"
        if sf.exists() and pf.exists():
            try:
                status = sf.read_text().strip()
                # Auto-fix stuck "running" jobs on list refresh
                if status == "running":
                    lf = d / "run.log"
                    if lf.exists():
                        tail = lf.read_text(errors="replace")[-500:]
                        if "Assessment complete" in tail:
                            sf.write_text("done"); status = "done"
                        elif "failed. (See above for error)" in tail:
                            sf.write_text("error"); status = "error"
                p = json.loads(pf.read_text())
                lf = d / "run.log"
                elapsed = round(time.time() - lf.stat().st_ctime) \
                          if lf.exists() else 0
                jobs.append({
                    "job_id":    d.name,
                    "status":    status,
                    "elapsed_s": elapsed,
                    **p
                })
            except Exception:
                pass
    return jobs



# ── Serve raw job files (PDB, PDBQT, SDF) for NGL viewer ─────────────
from fastapi.responses import FileResponse
import mimetypes

@app.get("/api/jobs/{job_id}/files/{file_path:path}")
def serve_job_file(job_id: str, file_path: str):
    full = JOBS_DIR / job_id / file_path
    if not full.exists() or not full.is_file():
        raise HTTPException(404, f"File not found: {file_path}")
    # Security: only allow known safe extensions
    allowed = {".pdb",".pdbqt",".sdf",".mol2",".log",".txt",".json",".pqr"}
    if full.suffix.lower() not in allowed:
        raise HTTPException(403, "File type not allowed")
    mt, _ = mimetypes.guess_type(str(full))
    return FileResponse(str(full), media_type=mt or "text/plain")


# ── 2D LigPlot-style interaction SVG ─────────────────────────────────
@app.get("/api/jobs/{job_id}/ligplot", response_class=HTMLResponse)
def job_ligplot(job_id: str):
    """Generate a 2D LigPlot-style SVG using RDKit + interactions.json"""
    jdir  = JOBS_DIR / job_id / "results"
    inter_f = jdir / "interactions" / "interactions.json"
    lig_f   = jdir / "ligand" / "compound_3D.sdf"
    ml_f    = jdir / "ml" / "ml_dtss_results.json"
    if not inter_f.exists():
        raise HTTPException(404, "interactions.json not found")

    try:
        from rdkit import Chem
        from rdkit.Chem import Draw, rdDepictor, AllChem
        from rdkit.Chem.Draw import rdMolDraw2D
        import json, math

        inter  = json.loads(inter_f.read_text())
        inters = inter.get("interactions", [])
        hbonds = [x for x in inters if x.get("type","").lower() in
                  ("hbond","h-bond","hydrogen bond","hbond_donor","hbond_acceptor")]
        hydro  = [x for x in inters if "hydro" in x.get("type","").lower() or
                  "hydrophobic" in x.get("type","").lower()]

        # Load or build ligand mol
        mol = None
        if lig_f.exists():
            mol = Chem.SDMolSupplier(str(lig_f), removeHs=True)[0]
        if mol is None:
            # Fallback: read SMILES from params
            params_f = JOBS_DIR / job_id / "params.json"
            if params_f.exists():
                params = json.loads(params_f.read_text())
                smi = params.get("smiles","")
                if smi:
                    mol = Chem.MolFromSmiles(smi)

        if mol is None:
            raise HTTPException(500, "Cannot load ligand structure")

        # Generate 2D coords
        rdDepictor.Compute2DCoords(mol)

        # ── Build SVG manually (LigPlot style) ─────────────────────
        # Get 2D coords of ligand atoms
        conf = mol.GetConformer()
        xs = [conf.GetAtomPosition(i).x for i in range(mol.GetNumAtoms())]
        ys = [conf.GetAtomPosition(i).y for i in range(mol.GetNumAtoms())]

        # Normalise to SVG coordinates (centre of canvas at 300,300)
        cx, cy = sum(xs)/len(xs), sum(ys)/len(ys)
        scale = 40
        SVG_CX, SVG_CY = 320, 320
        def to_svg(x, y):
            return (SVG_CX + (x - cx)*scale, SVG_CY - (y - cy)*scale)

        W, H = 700, 640
        svg_parts = [
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" '
            f'viewBox="0 0 {W} {H}">',
            '<rect width="100%" height="100%" fill="#fafbfc"/>',
            # Title
            '<text x="14" y="22" font-family="Arial,sans-serif" font-size="13" '
            'font-weight="bold" fill="#1a2a4a">Ligand–Protein Interaction Map</text>',
            '<text x="14" y="38" font-family="Arial,sans-serif" font-size="10" '
            'fill="#6b7c93">H-bonds: blue dashed  ·  Hydrophobic: green  ·  '
            'Best docking pose</text>',
        ]

        # ── Draw bonds ───────────────────────────────────────────
        for bond in mol.GetBonds():
            i1, i2 = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            x1,y1 = to_svg(xs[i1], ys[i1])
            x2,y2 = to_svg(xs[i2], ys[i2])
            btype = bond.GetBondTypeAsDouble()
            if btype >= 2:
                # Double bond — draw two lines
                dx, dy = y2-y1, x1-x2
                dl = max(math.sqrt(dx*dx+dy*dy), 0.01)
                ox, oy = dx/dl*2.5, dy/dl*2.5
                svg_parts.append(
                    f'<line x1="{x1+ox:.1f}" y1="{y1+oy:.1f}" '
                    f'x2="{x2+ox:.1f}" y2="{y2+oy:.1f}" '
                    'stroke="#444" stroke-width="1.4"/>')
                svg_parts.append(
                    f'<line x1="{x1-ox:.1f}" y1="{y1-oy:.1f}" '
                    f'x2="{x2-ox:.1f}" y2="{y2-oy:.1f}" '
                    'stroke="#444" stroke-width="1.4"/>')
            else:
                svg_parts.append(
                    f'<line x1="{x1:.1f}" y1="{y1:.1f}" '
                    f'x2="{x2:.1f}" y2="{y2:.1f}" '
                    'stroke="#444" stroke-width="1.8"/>')

        # ── Draw atoms ───────────────────────────────────────────
        ATOM_COLORS = {"O":"#e53935","N":"#1565c0","S":"#f9a825",
                       "P":"#e65100","F":"#2e7d32","Cl":"#558b2f",
                       "Br":"#6a1b9a","C":"#333","H":"#999"}
        for i in range(mol.GetNumAtoms()):
            atom  = mol.GetAtomWithIdx(i)
            sym   = atom.GetSymbol()
            ax,ay = to_svg(xs[i], ys[i])
            col   = ATOM_COLORS.get(sym, "#555")
            if sym != "C":
                svg_parts.append(
                    f'<circle cx="{ax:.1f}" cy="{ay:.1f}" r="10" '
                    f'fill="white" stroke="{col}" stroke-width="1.2"/>')
                svg_parts.append(
                    f'<text x="{ax:.1f}" y="{ay+4:.1f}" text-anchor="middle" '
                    f'font-family="Arial,sans-serif" font-size="10" '
                    f'font-weight="bold" fill="{col}">{sym}</text>')

        # ── Draw residue bubbles + interaction lines ──────────────
        # Arrange residues in a ring around the ligand
        seen = {}
        all_residues = []
        for ix in inters:
            res = ix.get("residue","")
            if res and res not in seen:
                seen[res] = ix.get("type","Hydrophobic")
                all_residues.append(res)

        n_res = len(all_residues)
        RING_R = 230
        for idx, res in enumerate(all_residues):
            angle  = (2*math.pi * idx / max(n_res,1)) - math.pi/2
            rx     = SVG_CX + RING_R * math.cos(angle)
            ry     = SVG_CY + RING_R * math.sin(angle)
            rtype  = seen[res]
            is_hb  = any(r.get("residue")==res for r in hbonds)
            is_hy  = any(r.get("residue")==res for r in hydro)
            rbg    = "#dbeafe" if is_hb else "#dcfce7"
            rbord  = "#1d4ed8" if is_hb else "#16a34a"
            rcol   = "#1d4ed8" if is_hb else "#15803d"
            lstyle = 'stroke-dasharray="6,3"' if is_hb else ''
            lstroke= "#1d4ed8" if is_hb else "#16a34a"
            lw     = "1.8" if is_hb else "1.4"

            # Interaction line to closest ligand atom
            lx, ly = SVG_CX, SVG_CY
            dist_f = next((r.get("distance",3.5) for r in inters
                          if r.get("residue")==res), 3.5)

            svg_parts.append(
                f'<line x1="{lx:.1f}" y1="{ly:.1f}" '
                f'x2="{rx:.1f}" y2="{ry:.1f}" '
                f'stroke="{lstroke}" stroke-width="{lw}" '
                f'opacity="0.6" {lstyle}/>')

            # Distance label
            mx, my = (lx+rx)/2, (ly+ry)/2
            svg_parts.append(
                f'<text x="{mx:.1f}" y="{my:.1f}" text-anchor="middle" '
                f'font-family="Arial,sans-serif" font-size="9" fill="{rcol}" '
                f'opacity="0.85">{dist_f:.2f}Å</text>')

            # Residue bubble
            svg_parts.append(
                f'<rect x="{rx-28:.1f}" y="{ry-12:.1f}" width="56" height="22" '
                f'rx="11" fill="{rbg}" stroke="{rbord}" stroke-width="1.5"/>')
            svg_parts.append(
                f'<text x="{rx:.1f}" y="{ry+5:.1f}" text-anchor="middle" '
                f'font-family="Arial,sans-serif" font-size="10" '
                f'font-weight="600" fill="{rcol}">{res}</text>')

        # ── Legend ───────────────────────────────────────────────
        lx0 = 14
        svg_parts += [
            f'<rect x="{lx0}" y="{H-52}" width="12" height="12" rx="3" '
            'fill="#dbeafe" stroke="#1d4ed8" stroke-width="1.5"/>',
            f'<text x="{lx0+17}" y="{H-42}" font-family="Arial,sans-serif" '
            f'font-size="11" fill="#444">H-bond ({len(hbonds)})</text>',
            f'<rect x="{lx0+110}" y="{H-52}" width="12" height="12" rx="3" '
            'fill="#dcfce7" stroke="#16a34a" stroke-width="1.5"/>',
            f'<text x="{lx0+127}" y="{H-42}" font-family="Arial,sans-serif" '
            f'font-size="11" fill="#444">Hydrophobic ({len(hydro)})</text>',
            f'<text x="{lx0}" y="{H-22}" font-family="Arial,sans-serif" '
            f'font-size="10" fill="#888">'
            f'Total interactions: {len(inters)}  ·  Residues: {n_res}</text>',
        ]

        svg_parts.append('</svg>')
        svg = '\n'.join(svg_parts)

        from fastapi.responses import HTMLResponse
        return HTMLResponse(content=svg, media_type="image/svg+xml")

    except ImportError as e:
        raise HTTPException(500, f"RDKit not available: {e}")
    except Exception as e:
        import traceback
        raise HTTPException(500, f"LigPlot error: {e}\n{traceback.format_exc()}")


# ── 2D Interaction Diagram — PNG (publication quality) ───────────────
@app.get("/api/jobs/{job_id}/interaction_png")
def interaction_png(job_id: str, download: int = 0):
    """Server-side RDKit 2D ligand-protein interaction PNG — pub-quality"""
    jdir   = JOBS_DIR / job_id / "results"
    # Return cached PNG if it exists
    cached = jdir / "interactions" / "interaction_2d.png"
    if cached.exists():
        from fastapi.responses import FileResponse
        headers = {}
        if download:
            headers["Content-Disposition"] = "attachment; filename=interaction_2d.png"
        return FileResponse(str(cached),
                            media_type="image/png", headers=headers)
    # Generate it now
    return _generate_interaction_png(job_id, jdir, cached, download=download)


def _generate_interaction_png(job_id, jdir, out_path, download=0):
    import json, math, io
    from pathlib import Path as P

    inter_f  = jdir / "interactions" / "interactions.json"
    params_f = JOBS_DIR / job_id / "params.json"
    ml_f     = jdir / "ml" / "ml_dtss_results.json"

    if not inter_f.exists():
        raise HTTPException(404, "interactions.json not found")

    try:
        from PIL import Image, ImageDraw, ImageFont
    except ImportError:
        raise HTTPException(500, "Pillow not installed: pip install Pillow")

    inter   = json.loads(inter_f.read_text())
    inters  = inter.get("interactions", [])

    hbond_res = {x["residue"] for x in inters
                 if any(k in x.get("type","").lower()
                        for k in ("hbond","h-bond","hydrogen"))}
    hydro_res = {x["residue"] for x in inters
                 if "hydro" in x.get("type","").lower()}
    all_res   = list({x["residue"]: x for x in inters}.keys())  # dedup ordered
    dist_map  = {x["residue"]: x.get("distance", 3.5)
                 for x in inters if x.get("residue")}

    # SMILES for ligand 2D structure via RDKit
    lig_svg_str = None
    try:
        from rdkit import Chem
        from rdkit.Chem import rdDepictor
        from rdkit.Chem.Draw import rdMolDraw2D

        mol = None
        lig_sdf = jdir / "ligand" / "compound_3D.sdf"
        if lig_sdf.exists():
            suppl = Chem.SDMolSupplier(str(lig_sdf), removeHs=True)
            mol   = suppl[0] if suppl else None
        if mol is None and params_f.exists():
            smi = json.loads(params_f.read_text()).get("smiles","")
            if smi: mol = Chem.MolFromSmiles(smi)
        if mol:
            rdDepictor.Compute2DCoords(mol)
            d = rdMolDraw2D.MolDraw2DSVG(300, 260)
            d.drawOptions().bondLineWidth = 2.5
            d.DrawMolecule(mol)
            d.FinishDrawing()
            lig_svg_str = d.GetDrawingText()
    except Exception:
        pass

    # Canvas size
    W, H   = 1000, 820
    img    = Image.new("RGB", (W, H), (250, 251, 252))
    draw   = ImageDraw.Draw(img, "RGBA")

    # Load fonts (use default PIL if not available)
    try:
        font_b = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 15)
        font_r = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 13)
        font_s = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 11)
    except Exception:
        font_b = ImageFont.load_default()
        font_r = font_b
        font_s = font_b

    # ── Draw ligand structure (SVG→PIL via cairosvg or placeholder) ──
    LIG_X, LIG_Y, LIG_W, LIG_H = 350, 260, 300, 260
    lig_placed = False

    if lig_svg_str:
        try:
            import cairosvg
            png_bytes = cairosvg.svg2png(bytestring=lig_svg_str.encode(),
                                         output_width=LIG_W, output_height=LIG_H)
            lig_img = Image.open(io.BytesIO(png_bytes)).convert("RGBA")
            # White background for ligand area
            bg = Image.new("RGBA", lig_img.size, (255,255,255,255))
            bg.paste(lig_img, mask=lig_img)
            img.paste(bg.convert("RGB"), (LIG_X, LIG_Y))
            draw.rounded_rectangle([LIG_X-2, LIG_Y-2, LIG_X+LIG_W+2, LIG_Y+LIG_H+2],
                                    radius=8, outline=(200,210,220,255), width=2)
            lig_placed = True
        except Exception:
            pass

    if not lig_placed:
        # Draw molecule name box as placeholder
        draw.rounded_rectangle([LIG_X, LIG_Y, LIG_X+LIG_W, LIG_Y+LIG_H],
                                radius=8, fill=(255,255,255,255),
                                outline=(26,42,74,255), width=2)
        draw.text((LIG_X+LIG_W//2, LIG_Y+LIG_H//2-10), "LIGAND",
                  fill=(26,42,74), font=font_b, anchor="mm")
        draw.text((LIG_X+LIG_W//2, LIG_Y+LIG_H//2+12), "(structure)",
                  fill=(150,160,170), font=font_s, anchor="mm")

    # ── Draw interaction ring of residues ──────────────────────────────
    cx, cy  = LIG_X + LIG_W//2, LIG_Y + LIG_H//2
    ring_r  = 310
    n_res   = len(all_res)

    for idx, res in enumerate(all_res):
        if not res: continue
        angle = (2 * math.pi * idx / max(n_res, 1)) - math.pi/2
        rx    = int(cx + ring_r * math.cos(angle))
        ry    = int(cy + ring_r * math.sin(angle))
        is_hb = res in hbond_res
        dist  = dist_map.get(res, 3.5)

        # Colours
        if is_hb:
            line_col  = (29, 78, 216, 200)
            fill_col  = (219, 234, 254, 230)
            text_col  = (29, 78, 216)
            bord_col  = (29, 78, 216, 255)
        else:
            line_col  = (22, 163, 74, 180)
            fill_col  = (220, 252, 231, 220)
            text_col  = (15, 100, 50)
            bord_col  = (22, 163, 74, 255)

        # Connection line
        if is_hb:
            # Dashed
            dx, dy = rx - cx, ry - cy
            dl     = max(math.sqrt(dx*dx + dy*dy), 1)
            ux, uy = dx/dl, dy/dl
            t = 0
            on = True
            while t < dl:
                if on:
                    x1 = int(cx + ux*t);     y1 = int(cy + uy*t)
                    x2 = int(cx + ux*min(t+9,dl)); y2 = int(cy + uy*min(t+9,dl))
                    draw.line([(x1,y1),(x2,y2)], fill=line_col, width=2)
                t += 12; on = not on
        else:
            draw.line([(cx,cy),(rx,ry)], fill=line_col, width=2)

        # Distance label midpoint
        mx, my = (cx+rx)//2, (cy+ry)//2
        draw.text((mx, my-7), f"{dist:.1f}Å",
                  fill=(100,110,120), font=font_s, anchor="mm")

        # Residue bubble
        bw, bh = 64, 26
        draw.rounded_rectangle([rx-bw//2, ry-bh//2, rx+bw//2, ry+bh//2],
                                radius=13, fill=fill_col, outline=bord_col, width=2)
        draw.text((rx, ry+1), res, fill=text_col, font=font_r, anchor="mm")

    # ── Title + stats ─────────────────────────────────────────────────
    draw.text((14, 14), "SAIDock — Ligand–Protein Interaction Map",
              fill=(26,42,74), font=font_b)

    ml_line = f"H-bonds: {len(hbond_res)}   Hydrophobic: {len(hydro_res)}   Total residues: {n_res}"
    if ml_f.exists():
        try:
            ml = json.loads(ml_f.read_text())
            ml_line += (f"   |   DTSS={ml.get('DTSS',0):.3f}"
                        f"   ΔG={ml.get('best_docking_score',0):.3f} kcal/mol")
        except Exception:
            pass
    draw.text((14, 36), ml_line, fill=(90,105,125), font=font_s)

    # ── Legend ────────────────────────────────────────────────────────
    ly = H - 44
    draw.rounded_rectangle([14, ly, 30, ly+16], radius=3,
                            fill=(219,234,254,255), outline=(29,78,216,255), width=2)
    draw.text((38, ly+2), f"H-bond ({len(hbond_res)})", fill=(29,78,216), font=font_s)
    draw.rounded_rectangle([160, ly, 176, ly+16], radius=3,
                            fill=(220,252,231,255), outline=(22,163,74,255), width=2)
    draw.text((184, ly+2), f"Hydrophobic ({len(hydro_res)})", fill=(22,163,74), font=font_s)
    draw.text((14, H-20), "Generated by SAIDock v1.0 | SSSIHL",
              fill=(180,190,200), font=font_s)

    # ── Save 300 DPI PNG ──────────────────────────────────────────────
    out_path.parent.mkdir(parents=True, exist_ok=True)
    buf = io.BytesIO()
    img.save(buf, format="PNG", dpi=(300, 300))
    buf.seek(0)
    out_path.write_bytes(buf.getvalue())

    from fastapi.responses import Response
    headers = {}
    if download:
        headers["Content-Disposition"] = "attachment; filename=interaction_2d.png"
    return Response(content=buf.getvalue(), media_type="image/png", headers=headers)


# ── Interaction SVG endpoint ──────────────────────────────────────────
@app.get("/api/jobs/{job_id}/interaction_svg")
def interaction_svg(job_id: str):
    from fastapi.responses import Response as FR
    jdir   = JOBS_DIR / job_id / "results"
    cached = jdir / "interactions" / "interaction_2d.svg"
    if cached.exists():
        return FR(content=cached.read_bytes(), media_type="image/svg+xml",
                  headers={"Content-Disposition": "attachment; filename=interaction_2d.svg"})
    raise HTTPException(404, "SVG not yet generated — view PNG first")


# ── Best pose PDBQT download ──────────────────────────────────────────
@app.get("/api/jobs/{job_id}/best_pose_pdbqt")
def best_pose_pdbqt(job_id: str):
    jdir  = JOBS_DIR / job_id / "results"
    ml_f  = jdir / "ml" / "ml_dtss_results.json"
    # Find best pocket id
    pocket_id = 1
    try:
        import json
        ml = json.loads(ml_f.read_text())
        # docking_results is in saidock_state.json
        state_f = jdir / "saidock_state.json"
        if state_f.exists():
            state = json.loads(state_f.read_text())
            dr = state.get("docking_results", [])
            if dr:
                best = min(dr, key=lambda x: x.get("best_score", 0))
                pocket_id = best.get("pocket_id", 1)
    except Exception:
        pass
    pdbqt = jdir / "docking" / f"pocket{pocket_id}_poses.pdbqt"
    if not pdbqt.exists():
        # Try any pdbqt
        pdbqts = list((jdir/"docking").glob("*_poses.pdbqt"))
        if pdbqts:
            pdbqt = sorted(pdbqts)[0]
        else:
            raise HTTPException(404, "No PDBQT pose found")
    from fastapi.responses import FileResponse
    return FileResponse(str(pdbqt), media_type="chemical/x-pdbqt",
                        filename=f"best_pose_pocket{pocket_id}.pdbqt",
                        headers={"Content-Disposition":
                                 f"attachment; filename=best_pose_pocket{pocket_id}.pdbqt"})


# ── Serve static frontend ─────────────────────────────────────────────
_STATIC = Path("web/static")
_STATIC.mkdir(parents=True, exist_ok=True)

try:
    app.mount("/static", StaticFiles(directory=str(_STATIC)), name="static")
except Exception:
    pass

@app.get("/", response_class=HTMLResponse)
def root():
    idx = _STATIC / "index.html"
    return HTMLResponse(idx.read_text() if idx.exists() else
                        "<h2>SAIDock: index.html not found</h2>")

@app.get("/api/jobs/{job_id}/report", response_class=HTMLResponse)
def job_report(job_id: str):
    """Serve HTML report with no-cache headers."""
    report = _find_report(JOBS_DIR / job_id)
    if not report:
        raise HTTPException(status_code=404, detail="Report not found — run a job first")
    return HTMLResponse(
        content=Path(report).read_text(encoding="utf-8"),
        headers={
            "Cache-Control": "no-store, no-cache, must-revalidate, max-age=0",
            "Pragma": "no-cache",
            "Expires": "0",
        }
    )
