# -*- coding: utf-8 -*-
import sys
if sys.version_info < (3, 6):
    sys.exit(
        "ERROR: Python 3.6+ required. Your version: {}\n\n"
        "On Windows run:  py -3 vina_docking.py ...".format(sys.version)
    )
"""
Vina Batch Docking — Parallel Ligand Docking Automation
=========================================================
Runs AutoDock Vina for all ligands in a folder against a single receptor.
Saves the best binding affinity for each ligand to a CSV file.

Features:
  - Parallel execution (up to 10 ligands simultaneously)
  - Real-time progress display per ligand
  - Best affinity extracted and saved to CSV
  - Graceful error handling per ligand

Usage:
    python vina_docking.py --ligands  ./ligands --receptor ./receptor.pdbqt --output   ./affinity.csv --config   ./config.txt --workers  10

Config file (config.txt) example:
    center_x = -21.0
    center_y =  -0.9
    center_z = -47.4
    size_x   =  20
    size_y   =  20
    size_z   =  20
    exhaustiveness = 8
    num_modes = 9
    energy_range = 3

Dependencies:
    AutoDock Vina must be installed and in PATH.
    pip install rich   (for the progress display)
"""

import argparse
import csv
import logging
import re
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path
from threading import Lock

# ── Optional rich for pretty progress ──────────────────────────────────────
try:
    from rich.console import Console
    from rich.live import Live
    from rich.panel import Panel
    from rich.progress import (
        BarColumn,
        MofNCompleteColumn,
        Progress,
        SpinnerColumn,
        TaskID,
        TextColumn,
        TimeElapsedColumn,
    )
    from rich.table import Table
    from rich.text import Text
    RICH_AVAILABLE = True
except ImportError:
    RICH_AVAILABLE = False

logging.basicConfig(level=logging.INFO, format="%(message)s")
log = logging.getLogger(__name__)

SUPPORTED_EXT = {".pdbqt"}

# ── CSV lock (thread-safe writes) ───────────────────────────────────────────
csv_lock = Lock()


# ============================================================
#  Config helpers
# ============================================================

def parse_config(config_path):
    """Parse a Vina config file into a dict of key->value strings."""
    config = {}
    if config_path is None:
        return config
    path = Path(config_path)
    if not path.exists():
        log.error("Config file not found: %s", config_path)
        sys.exit(1)
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if "=" in line:
                key, _, val = line.partition("=")
                config[key.strip()] = val.strip()
    return config


def build_vina_command(vina_bin, receptor, ligand, output_path, config, extra_args):
    """Assemble the vina command list."""
    cmd = [
        vina_bin,
        "--receptor", str(receptor),
        "--ligand",   str(ligand),
        "--out",      str(output_path),
    ]
    # Append config values as CLI flags
    for key, val in config.items():
        cmd += ["--{}".format(key), val]
    # Extra pass-through args
    if extra_args:
        cmd += extra_args
    return cmd


def extract_best_affinity(vina_output):
    """
    Parse Vina output and return the best affinity value (kcal/mol).
    Handles Vina 1.1.x and 1.2.x output formats.

    Vina 1.2.x sometimes outputs 0 for mode 1 even on success —
    the real best score appears in subsequent modes. This function
    scans all modes and returns the first non-zero value found.
    If ALL modes are 0, returns None (docking truly failed).
    """
    lines = vina_output.splitlines()
    in_table = False
    best = None
    for line in lines:
        if re.search(r"-{3,}", line):
            in_table = True
            continue
        if not in_table:
            continue
        parts = line.replace("|", " ").split()
        if not parts:
            continue
        try:
            int(parts[0])  # first token must be mode number
        except ValueError:
            continue
        if len(parts) < 2:
            continue
        try:
            val = float(parts[1])
            if val != 0.0:
                # Return immediately — modes are sorted best-first,
                # so first non-zero is the best valid score
                return val
            # val == 0: skip this mode, check the next ones
        except (ValueError, IndexError):
            continue
    # All modes were 0 — docking failed silently
    return None


# ============================================================
#  Ligand pre-processing (fix common format issues)
# ============================================================

def _has_valid_charges(atom_lines):
    """
    Check if ATOM/HETATM lines have non-zero partial charges.
    Charge is the second-to-last field in a PDBQT line.
    Returns True if at least half the atoms have non-zero charges.
    """
    nonzero = 0
    total = 0
    for line in atom_lines:
        parts = line.split()
        if len(parts) < 2:
            continue
        try:
            charge = float(parts[-2])
            if charge != 0.0:
                nonzero += 1
            total += 1
        except (ValueError, IndexError):
            continue
    return total > 0 and (nonzero / total) > 0.5


def _obabel_regenerate(ligand_path, tmp_dir):
    """
    Use Open Babel to regenerate a proper PDBQT with:
      - Gasteiger partial charges
      - Rotatable bond detection (BRANCH/TORSDOF)
    Returns path to new file, or None if obabel is not available.
    """
    import tempfile, os
    tmp_fd, tmp_path = tempfile.mkstemp(
        suffix=".pdbqt",
        prefix="vinafix_{}_".format(ligand_path.stem),
        dir=tmp_dir,
    )
    os.close(tmp_fd)
    out = Path(tmp_path)
    try:
        result = subprocess.run(
            [
                "obabel",
                "-ipdbqt", str(ligand_path),
                "-opdbqt",
                "-O", str(out),
                "--partialcharge", "gasteiger",
            ],
            capture_output=True, text=True, timeout=60,
        )
        if result.returncode == 0 and out.exists() and out.stat().st_size > 0:
            return out
        out.unlink(missing_ok=True)
        return None
    except (FileNotFoundError, subprocess.TimeoutExpired):
        out.unlink(missing_ok=True)
        return None


def preprocess_ligand(ligand_path, tmp_dir):
    """
    Inspect a ligand PDBQT and fix issues that cause Vina to return
    score=0 or reject the file outright.

    Checks performed (in order):
      1. Double CRLF (\r\r\n) — normalize to LF.
      2. Invalid tags (CONECT, MASTER, …) — strip them.
      3. Missing charges — if all/most charges are 0, attempt to
         recalculate with Open Babel (Gasteiger).
      4. Missing rotatable bonds (TORSDOF 0 / no BRANCH) — attempt
         to regenerate torsion tree with Open Babel.
         If obabel is unavailable, add a rigid wrapper as fallback
         (score will be 0 but docking won't crash).

    Returns the path to use for docking (original or fixed temp file).
    Caller is responsible for deleting any temp file created.
    """
    import tempfile, os

    raw = ligand_path.read_bytes()

    # ── Fix 1: double CRLF ──────────────────────────────────
    fixed = raw.replace(b"\r\r\n", b"\n").replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    text = fixed.decode("utf-8", errors="replace")

    # ── Fix 2: strip tags Vina does not accept in ligands ───
    ALLOWED = ("REMARK", "ATOM", "HETATM", "ROOT", "ENDROOT",
               "BRANCH", "ENDBRANCH", "TORSDOF", "MODEL", "ENDMDL", "TER")
    lines = text.splitlines()
    filtered = [l for l in lines if not l.strip() or l.startswith(ALLOWED)]
    removed = len(lines) - len(filtered)
    lines = filtered
    text = "\n".join(lines) + "\n"

    atom_lines = [l for l in lines if l.startswith(("ATOM", "HETATM"))]
    has_root    = any(l.startswith("ROOT") for l in lines)
    has_torsdof = any(l.startswith("TORSDOF") for l in lines)
    has_branch  = any(l.startswith("BRANCH") for l in lines)
    torsdof_val = 0
    for l in lines:
        if l.startswith("TORSDOF"):
            try:
                torsdof_val = int(l.split()[1])
            except (IndexError, ValueError):
                pass

    needs_obabel = (
        not _has_valid_charges(atom_lines)   # charges missing/zero
        or not has_branch                     # no rotatable bonds defined
        or not has_root                       # no torsion tree at all
        or torsdof_val == 0                   # explicitly rigid
    )

    if needs_obabel:
        # Write current cleaned version to temp so obabel can read it
        tmp_fd, clean_path = tempfile.mkstemp(
            suffix=".pdbqt", prefix="vinaclean_{}_".format(ligand_path.stem), dir=tmp_dir
        )
        with os.fdopen(tmp_fd, "w", encoding="utf-8") as fh:
            fh.write(text)
        clean = Path(clean_path)

        regen = _obabel_regenerate(clean, tmp_dir)
        clean.unlink(missing_ok=True)

        if regen is not None:
            return regen  # obabel produced a proper flexible PDBQT

        # obabel not available — fall back to rigid wrapper
        log.warning(
            "  [WARN] %s: no rotatable bonds / charges detected and obabel is "
            "not available. Docking as rigid molecule (score may be 0). "
            "Install obabel for proper results.", ligand_path.name
        )
        remark_lines = [l for l in lines if l.startswith("REMARK")]
        rebuilt = remark_lines + ["ROOT"] + atom_lines + ["ENDROOT", "TORSDOF 0"]
        text = "\n".join(rebuilt) + "\n"

    elif removed > 0 or fixed != raw:
        pass  # text already updated above
    else:
        return ligand_path  # nothing changed — use original

    tmp_fd, tmp_path = tempfile.mkstemp(
        suffix=".pdbqt",
        prefix="vinafix_{}_".format(ligand_path.stem),
        dir=tmp_dir,
    )
    with os.fdopen(tmp_fd, "w", encoding="utf-8") as fh:
        fh.write(text)
    return Path(tmp_path)


# ============================================================
#  Single docking job
# ============================================================

def run_docking(
    ligand_path,
    receptor_path,
    vina_bin,
    config,
    extra_args,
    progress=None,
    task_id=None,
    status_dict=None,
    status_lock=None,
    debug=False,
):
    """
    Run a single Vina docking job.
    Returns (ligand_name, affinity_kcal, elapsed_s, error_msg).
    Pose output is written to a temp file and deleted immediately after.
    """
    import tempfile, os
    name = ligand_path.stem
    tmp_dir = Path(tempfile.gettempdir())

    # Pre-process ligand: fix CRLF and missing torsion tree if needed
    fixed_ligand = preprocess_ligand(ligand_path, tmp_dir)
    ligand_was_fixed = fixed_ligand != ligand_path

    tmp_fd, tmp_path = tempfile.mkstemp(suffix=".pdbqt", prefix="vina_{}_".format(name))
    os.close(tmp_fd)
    out_file = Path(tmp_path)
    cmd = build_vina_command(
        vina_bin, receptor_path, fixed_ligand, out_file, config, extra_args
    )

    def _set_status(msg):
        if status_dict is not None and status_lock is not None:
            with status_lock:
                status_dict[name] = msg
        if progress is not None and task_id is not None:
            progress.update(task_id, description=f"[cyan]{name:<30}[/cyan] {msg}")

    _set_status("starting…")
    t0 = time.monotonic()

    try:
        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,  # merge stderr into stdout — Vina writes results to stderr
            text=True,
            timeout=3600,
        )
        elapsed = time.monotonic() - t0

        # Vina 1.2.x sometimes exits with non-zero even on success.
        # So we always try to parse affinity first, and only treat
        # it as a real failure if parsing also comes up empty.
        affinity = extract_best_affinity(result.stdout)

        if affinity is None:
            if result.returncode != 0:
                lines = result.stdout.strip().splitlines()
                error_lines = [
                    l for l in lines
                    if l.strip()
                    and not l.startswith("#")
                    and "AutoDock Vina" not in l
                    and "cite" not in l.lower()
                    and "http" not in l
                    and l.strip() != ""
                ]
                err = " | ".join(error_lines[:5]) if error_lines else result.stdout.strip()[-300:]
            else:
                err = "Could not parse affinity from Vina output."
            if debug:
                print("\n[DEBUG] Full Vina output for {}:\n{}".format(name, result.stdout))
            _set_status(f"[red]FAILED[/red]")
            return name, None, elapsed, err, result.stdout

        _set_status(
            f"[green]done[/green] [bold]{affinity:.2f} kcal/mol[/bold] "
            f"[dim]({elapsed:.0f}s)[/dim]"
        )
        if progress and task_id is not None:
            progress.advance(task_id)
        return name, affinity, elapsed, None, result.stdout

    except FileNotFoundError:
        elapsed = time.monotonic() - t0
        _set_status("[red]vina not found[/red]")
        return name, None, elapsed, (
            f"'{vina_bin}' not found. Install Vina and ensure it is in PATH."
        ), ""
    except subprocess.TimeoutExpired:
        elapsed = time.monotonic() - t0
        _set_status("[red]timeout[/red]")
        return name, None, elapsed, "Vina timed out after 3600 s.", ""
    except Exception as exc:
        elapsed = time.monotonic() - t0
        _set_status(f"[red]error: {exc}[/red]")
        return name, None, elapsed, str(exc), ""
    finally:
        # Delete temp pose file — we only need the affinity value
        try:
            out_file.unlink(missing_ok=True)
        except Exception:
            pass
        # Delete fixed ligand temp file if one was created
        if ligand_was_fixed:
            try:
                fixed_ligand.unlink(missing_ok=True)
            except Exception:
                pass


# ============================================================
#  Progress display (rich)
# ============================================================

def _make_progress():
    return Progress(
        SpinnerColumn(),
        TextColumn("{task.description}"),
        BarColumn(bar_width=20),
        MofNCompleteColumn(),
        TimeElapsedColumn(),
        expand=False,
    )


# ============================================================
#  CSV writer
# ============================================================

def append_csv_row(csv_path, row):
    with csv_lock:
        with open(csv_path, "a", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh)
            writer.writerow(row)


def init_csv(csv_path):
    with open(csv_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["ligand", "affinity_kcal_mol", "elapsed_s", "status"])


# ============================================================
#  Debug log writer
# ============================================================

def write_debug_log(debug_dir, name, affinity, elapsed, err, vina_output,
                    session_log_lock=None, session_log_path=None):
    """
    Write per-ligand debug info in two places:
      1. debug_logs/<name>.log  — individual file per ligand
      2. debug_logs/session.log — consolidated log for the whole run
    """
    debug_dir.mkdir(parents=True, exist_ok=True)
    status = "OK" if affinity is not None else "FAILED"

    def _block(fh):
        fh.write("=" * 60 + "\n")
        fh.write("  Ligand   : {}\n".format(name))
        fh.write("  Status   : {}\n".format(status))
        if affinity is not None:
            fh.write("  Affinity : {:.4g} kcal/mol\n".format(affinity))
        else:
            fh.write("  Error    : {}\n".format(err))
        fh.write("  Elapsed  : {:.1f} s\n".format(elapsed))
        fh.write("  Time     : {}\n".format(datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
        fh.write("=" * 60 + "\n\n")
        fh.write("── Vina Output ─────────────────────────────────────────\n")
        fh.write(vina_output if vina_output else "(no output captured)\n")
        fh.write("\n")

    # Individual log file
    with open(debug_dir / "{}.log".format(name), "w", encoding="utf-8") as fh:
        _block(fh)

    # Consolidated session log (thread-safe append)
    if session_log_path is not None:
        lock = session_log_lock if session_log_lock is not None else Lock()
        with lock:
            with open(session_log_path, "a", encoding="utf-8") as fh:
                _block(fh)




# ============================================================
#  Main orchestrator
# ============================================================

def resolve_ligands(inputs):
    """
    Accept any mix of:
      - a folder path            -> all .pdbqt inside it
      - individual .pdbqt files  -> used directly
      - a .txt file              -> each non-empty line treated as a path
    Returns a sorted, deduplicated list of Path objects.
    """
    seen = set()
    ligands = []
    for raw in inputs:
        p = Path(raw).resolve()
        # .txt list file
        if p.suffix.lower() == ".txt":
            if not p.exists():
                log.error("List file not found: %s", p)
                sys.exit(1)
            with open(p, encoding="utf-8", errors="ignore") as fh:
                for line in fh:
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue
                    fp = Path(line).resolve()
                    if not fp.exists():
                        log.warning("  [SKIP] Not found: %s", fp)
                        continue
                    if fp.suffix.lower() not in SUPPORTED_EXT:
                        log.warning("  [SKIP] Unsupported format: %s", fp.name)
                        continue
                    if fp not in seen:
                        seen.add(fp)
                        ligands.append(fp)
        # folder
        elif p.is_dir():
            for f in sorted(p.iterdir()):
                if f.is_file() and f.suffix.lower() in SUPPORTED_EXT:
                    if f not in seen:
                        seen.add(f)
                        ligands.append(f)
        # individual file
        elif p.is_file():
            if p.suffix.lower() not in SUPPORTED_EXT:
                log.warning("  [SKIP] Unsupported format: %s", p.name)
                continue
            if p not in seen:
                seen.add(p)
                ligands.append(p)
        else:
            log.warning("  [SKIP] Path not found: %s", p)
    return ligands


def run_batch(
    ligands_inputs,
    receptor,
    output_csv,
    vina_bin,
    config_path,
    workers,
    extra_args,
    debug=False,
):
    receptor_path = Path(receptor).resolve()
    csv_path = Path(output_csv).resolve()
    debug_dir = csv_path.parent / "debug_logs"

    # ── Validate inputs ─────────────────────────────────────
    if not receptor_path.exists():
        log.error("Receptor file not found: %s", receptor_path)
        sys.exit(1)

    ligands = resolve_ligands(ligands_inputs)
    if not ligands:
        log.error("No valid .pdbqt ligand files found in the provided inputs.")
        sys.exit(1)

    config = parse_config(config_path)
    init_csv(csv_path)

    # ── Debug log setup ─────────────────────────────────────
    debug_dir.mkdir(parents=True, exist_ok=True)
    session_log_path = debug_dir / "session.log"
    session_log_lock = Lock()
    # Write session header
    with open(session_log_path, "w", encoding="utf-8") as fh:
        fh.write("Vina Batch Docking — Session Log\n")
        fh.write("Started  : {}\n".format(datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
        fh.write("Receptor : {}\n".format(receptor_path.name))
        fh.write("Ligands  : {} files\n".format(len(ligands) if ligands else "?"))
        fh.write("=" * 60 + "\n\n")

    workers = min(workers, 10, len(ligands))

    # ── Header ──────────────────────────────────────────────
    sep = "═" * 60
    print(f"\n{sep}")
    print(f"  Vina Batch Docking")
    print(f"  {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"  Receptor : {receptor_path.name}")
    print(f"  Ligands  : {len(ligands)} files")
    print(f"  Output   : {csv_path}")
    print(f"  Debug    : {debug_dir}")
    print(f"  Workers  : {workers} parallel")
    if config:
        print(f"  Config   : {config_path}")
    print(f"{sep}\n")

    results = []
    status_dict = {}
    status_lock = Lock()

    # ── Rich progress display ────────────────────────────────
    if RICH_AVAILABLE:
        console = Console()
        progress = _make_progress()
        overall = progress.add_task(
            f"[bold white]Overall ({len(ligands)} ligands)", total=len(ligands)
        )
        tasks = {}
        for lig in ligands:
            tid = progress.add_task(
                f"[dim]{lig.stem:<30}[/dim] queued",
                total=1,
                visible=True,
            )
            tasks[lig.stem] = tid

        with Live(progress, console=console, refresh_per_second=8):
            with ThreadPoolExecutor(max_workers=workers) as executor:
                futures = {
                    executor.submit(
                        run_docking,
                        lig,
                        receptor_path,
                        vina_bin,
                        config,
                        extra_args,
                        progress,
                        tasks[lig.stem],
                        status_dict,
                        status_lock,
                        debug,
                    ): lig
                    for lig in ligands
                }
                for future in as_completed(futures):
                    name, affinity, elapsed, err, vina_out = future.result()
                    status = "ok" if affinity is not None else "error"
                    append_csv_row(csv_path, [name, affinity if affinity is not None else "", f"{elapsed:.1f}", status])
                    write_debug_log(debug_dir, name, affinity, elapsed, err, vina_out,
                                    session_log_lock, session_log_path)
                    results.append((name, affinity, elapsed, err))
                    progress.advance(overall)

    else:
        # ── Fallback: plain logging ──────────────────────────
        log.info("(Install 'rich' for a better progress display: pip install rich)\n")
        with ThreadPoolExecutor(max_workers=workers) as executor:
            futures = {
                executor.submit(
                    run_docking,
                    lig,
                    receptor_path,
                    vina_bin,
                    config,
                    extra_args,
                    None, None,
                    status_dict,
                    status_lock,
                    debug,
                ): lig
                for lig in ligands
            }
            for future in as_completed(futures):
                name, affinity, elapsed, err, vina_out = future.result()
                status = "ok" if affinity is not None else "error"
                append_csv_row(csv_path, [name, affinity if affinity is not None else "", f"{elapsed:.1f}", status])
                write_debug_log(debug_dir, name, affinity, elapsed, err, vina_out,
                                session_log_lock, session_log_path)
                results.append((name, affinity, elapsed, err))
                tag = f"{affinity:.2f} kcal/mol" if affinity else f"FAILED — {err}"
                log.info("  [%s] %s  %s  (%.0fs)", status.upper(), name, tag, elapsed)

    # ── Summary ──────────────────────────────────────────────
    ok = [(n, a, t) for n, a, t, e in results if a is not None]
    failed = [(n, e) for n, a, t, e in results if a is None]

    ok_sorted = sorted(ok, key=lambda x: x[1])

    print(f"\n{sep}")
    print(f"  Results saved → {csv_path}")
    print(f"  Completed: {len(ok)}/{len(ligands)}  |  Errors: {len(failed)}")
    if ok_sorted:
        best = ok_sorted[0]
        print(f"\n  🏆 Best affinity: {best[0]}  →  {best[1]:.2f} kcal/mol")
        print(f"\n  Top 5:")
        for i, (n, a, t) in enumerate(ok_sorted[:5], 1):
            print(f"    {i}. {n:<35} {a:>8.2f} kcal/mol")
    if failed:
        print(f"\n  Errors ({len(failed)}):")
        for n, e in failed:
            print(f"    ✗ {n}: {e}")
    print(f"{sep}\n")

    # Write session footer to consolidated log
    with open(session_log_path, "a", encoding="utf-8") as fh:
        fh.write("\n" + "=" * 60 + "\n")
        fh.write("Session Summary\n")
        fh.write("Finished : {}\n".format(datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
        fh.write("Completed: {}/{}\n".format(len(ok), len(results)))
        fh.write("Errors   : {}\n".format(len(failed)))
        if ok_sorted:
            fh.write("Best     : {}  {:.4g} kcal/mol\n".format(ok_sorted[0][0], ok_sorted[0][1]))
        fh.write("=" * 60 + "\n")


# ============================================================
#  CLI
# ============================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Parallel AutoDock Vina batch docking",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="""
Examples:
  # Folder of ligands (original usage):
  python vina_docking.py --ligands ./ligands --receptor receptor.pdbqt --config config.txt

  # Individual .pdbqt files:
  python vina_docking.py --ligands FDA_001.pdbqt NuBBE_1.pdbqt --receptor receptor.pdbqt --config config.txt

  # Mix of folder + individual files:
  python vina_docking.py --ligands ./ligands FDA_extra.pdbqt --receptor receptor.pdbqt --config config.txt

  # Text file with one path per line:
  python vina_docking.py --ligands my_list.txt --receptor receptor.pdbqt --config config.txt

  # Extra Vina args (after --):
  python vina_docking.py --ligands ./ligands --receptor rec.pdbqt --config config.txt -- --exhaustiveness 16
        """
    )

    parser.add_argument("--ligands", required=True, nargs="+",
                        help=(
                            "One or more inputs (space-separated):\n"
                            "  folder/      -> all .pdbqt inside\n"
                            "  file.pdbqt   -> individual ligand\n"
                            "  list.txt     -> text file with one path per line"
                        ))
    parser.add_argument("--receptor",  required=True,
                        help="Receptor .pdbqt file")
    parser.add_argument("--config",    default=None,
                        help="Vina config file (center_x/y/z, size_x/y/z, etc.)")
    parser.add_argument("--output",    default="affinity.csv",
                        help="Output CSV file (default: affinity.csv)")
    parser.add_argument("--workers",   type=int, default=4,
                        help="Parallel workers, max 10 (default: 4)")
    parser.add_argument("--vina",      default="vina",
                        help="Path/name of the Vina executable (default: vina)")
    parser.add_argument("--debug", action="store_true",
                        help="Print full Vina output for failed ligands")
    parser.add_argument("extra", nargs=argparse.REMAINDER,
                        help="Extra args passed directly to Vina (after --)")

    args = parser.parse_args()

    if args.workers > 10:
        log.warning("--workers capped at 10 (requested %d)", args.workers)
        args.workers = 10

    # Strip leading '--' separator if present
    extra = [a for a in (args.extra or []) if a != "--"]

    run_batch(
        ligands_inputs=args.ligands,
        receptor=args.receptor,
        output_csv=args.output,
        vina_bin=args.vina,
        config_path=args.config,
        workers=args.workers,
        extra_args=extra,
        debug=args.debug,
    )
