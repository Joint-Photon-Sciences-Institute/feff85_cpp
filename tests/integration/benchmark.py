#!/usr/bin/env python3
"""
Benchmark script for FEFF85 C++ vs Fortran baseline.

Runs feff8l.exe in each material's SCF and noSCF folders, compares
the resulting chi.dat against the Fortran baseline, and produces
summary comparison plots.
"""

import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt


# ---------- paths ----------------------------------------------------------

SCRIPT_DIR = Path(__file__).resolve().parent
WORK_DIR = SCRIPT_DIR / "work"
FEFF_EXE = SCRIPT_DIR.parent.parent / "build" / "apps" / "feff8l.exe"

MATERIALS = sorted(
    [d.name for d in WORK_DIR.iterdir() if d.is_dir()],
    key=str.lower,
)

# Mapping from run‐folder name to baseline subfolder name
BASELINE_MAP = {"SCF": "withSCF", "noSCF": "noSCF"}


# ---------- helpers --------------------------------------------------------

def parse_chi(path: Path) -> tuple[np.ndarray, np.ndarray]:
    """Read a chi.dat file and return (k, chi) arrays."""
    k_vals, chi_vals = [], []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 2:
                k_vals.append(float(parts[0]))
                chi_vals.append(float(parts[1]))
    return np.array(k_vals), np.array(chi_vals)


def r_factor(chi_test: np.ndarray, chi_ref: np.ndarray) -> float:
    """Compute R-factor between two chi arrays (must be same length)."""
    denom = np.sum(chi_ref ** 2)
    if denom == 0:
        return float("inf")
    return float(np.sum((chi_test - chi_ref) ** 2) / denom)


def run_feff(folder: Path) -> tuple[int, float]:
    """Run feff8l.exe inside *folder* and return (returncode, elapsed_s)."""
    t0 = time.perf_counter()
    result = subprocess.run(
        [str(FEFF_EXE)],
        cwd=str(folder),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    )
    elapsed = time.perf_counter() - t0
    if result.returncode != 0:
        print(f"  [FAIL] feff8l returned {result.returncode} in {folder}")
        print(result.stdout[-500:] if result.stdout else "")
    return result.returncode, elapsed


# ---------- main -----------------------------------------------------------

def main():
    if not FEFF_EXE.exists():
        sys.exit(f"ERROR: executable not found: {FEFF_EXE}")

    print(f"Executable : {FEFF_EXE}")
    print(f"Work dir   : {WORK_DIR}")
    print(f"Materials  : {', '.join(MATERIALS)}")
    print()

    # Collect results: results[mode][material] = (k_cpp, chi_cpp, k_ref, chi_ref, rfact, elapsed, rc)
    results: dict[str, dict] = {"SCF": {}, "noSCF": {}}

    for material in MATERIALS:
        for mode in ("noSCF", "SCF"):
            run_dir = WORK_DIR / material / mode
            baseline_dir = WORK_DIR / material / "baseline" / BASELINE_MAP[mode]

            if not run_dir.exists():
                print(f"  SKIP {material}/{mode} — run folder missing")
                continue
            if not (baseline_dir / "chi.dat").exists():
                print(f"  SKIP {material}/{mode} — baseline chi.dat missing")
                continue

            print(f"Running {material}/{mode} ... ", end="", flush=True)
            rc, elapsed = run_feff(run_dir)

            chi_path = run_dir / "chi.dat"
            if rc != 0 or not chi_path.exists():
                print(f"FAILED ({elapsed:.1f}s)")
                results[mode][material] = None
                continue

            k_cpp, chi_cpp = parse_chi(chi_path)
            k_ref, chi_ref = parse_chi(baseline_dir / "chi.dat")

            # Interpolate to common grid if lengths differ
            if len(k_cpp) != len(k_ref) or not np.allclose(k_cpp, k_ref, atol=1e-6):
                k_common = np.linspace(
                    max(k_cpp[0], k_ref[0]),
                    min(k_cpp[-1], k_ref[-1]),
                    min(len(k_cpp), len(k_ref)),
                )
                chi_cpp_i = np.interp(k_common, k_cpp, chi_cpp)
                chi_ref_i = np.interp(k_common, k_ref, chi_ref)
            else:
                k_common = k_cpp
                chi_cpp_i = chi_cpp
                chi_ref_i = chi_ref

            rf = r_factor(chi_cpp_i, chi_ref_i)
            print(f"done  R={rf:.2e}  ({elapsed:.1f}s)")

            results[mode][material] = {
                "k_cpp": k_common,
                "chi_cpp": chi_cpp_i,
                "k_ref": k_common,
                "chi_ref": chi_ref_i,
                "rfactor": rf,
                "elapsed": elapsed,
            }

    # ---------- plots ------------------------------------------------------

    for mode in ("noSCF", "SCF"):
        mode_label = "SCF" if mode == "SCF" else "noSCF"
        materials_with_data = [m for m in MATERIALS if results[mode].get(m)]

        if not materials_with_data:
            print(f"\nNo data for {mode_label} — skipping plot.")
            continue

        n = len(materials_with_data)
        ncols = 2
        nrows = (n + ncols - 1) // ncols
        fig, axes = plt.subplots(nrows, ncols, figsize=(14, 4 * nrows))
        axes = np.atleast_2d(axes)

        for idx, material in enumerate(materials_with_data):
            row, col = divmod(idx, ncols)
            ax = axes[row, col]
            d = results[mode][material]

            ax.plot(d["k_ref"], d["chi_ref"], label="Fortran baseline", linewidth=1.2)
            ax.plot(d["k_cpp"], d["chi_cpp"], label="C++ (feff8l)", linewidth=1.2,
                    linestyle="--")
            ax.set_title(f"{material}  (R = {d['rfactor']:.2e})")
            ax.set_xlabel("k (1/A)")
            ax.set_ylabel("chi(k)")
            ax.legend(fontsize=8)
            ax.grid(True, alpha=0.3)

        # Hide unused subplots
        for idx in range(n, nrows * ncols):
            row, col = divmod(idx, ncols)
            axes[row, col].set_visible(False)

        fig.suptitle(f"FEFF85 C++ vs Fortran — {mode_label}", fontsize=14, y=1.01)
        fig.tight_layout()

        out_path = SCRIPT_DIR / f"ALL_{mode_label}_final.png"
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"\nSaved: {out_path}")

    # ---------- summary table ----------------------------------------------

    print("\n" + "=" * 70)
    print(f"{'Material':<20} {'Mode':<8} {'R-factor':<14} {'Time (s)':<10} {'Status'}")
    print("-" * 70)
    for material in MATERIALS:
        for mode in ("noSCF", "SCF"):
            d = results[mode].get(material)
            if d is None:
                print(f"{material:<20} {mode:<8} {'—':<14} {'—':<10} FAIL")
            else:
                status = "PASS" if d["rfactor"] < 0.01 else "CHECK"
                print(f"{material:<20} {mode:<8} {d['rfactor']:<14.2e} {d['elapsed']:<10.1f} {status}")
    print("=" * 70)


if __name__ == "__main__":
    main()
