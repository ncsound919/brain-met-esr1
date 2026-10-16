"""Run the dormancy model across a parameter grid.

Pre-registered prediction: diversification risk rises under pause-targeting.
Falsification: flat or reversed effect across the plausible range.

Uses GillesPy2 SSA when installed; otherwise a NumPy tau-leap fallback so
smoke runs and tests need no extra dependency. Writes results/m4 CSV.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from model import build_model

try:
    from rich.console import Console
    from rich.progress import track
    from rich.table import Table

    HAVE_RICH = True
except ImportError:
    HAVE_RICH = False

ROOT = Path(__file__).resolve().parents[2]


def fallback_tauleap(
    switching_rate: float,
    pause_rate: float,
    seed: int = 0,
    t_end: int = 500,
    d0: int = 100,
) -> int:
    rng = np.random.default_rng(seed)
    D, A = d0, 0
    for _ in range(t_end):
        n_switch = min(rng.poisson(switching_rate * D), D)
        n_rev = min(rng.poisson(switching_rate * A), A)
        D = D - n_switch + n_rev
        A = A + n_switch - n_rev
        n_prolif = rng.poisson(0.1 * A)
        n_death = min(rng.poisson((0.05 + pause_rate) * A), A + n_prolif)
        A = A + n_prolif - n_death
        D, A = max(D, 0), max(A, 0)
    return int(A)


def run_grid(
    switching_rates: list[float],
    pause_rates: list[float],
    n_trajectories: int = 50,
    t_end: int = 500,
    show_progress: bool = True,
) -> pd.DataFrame:
    try:
        import gillespy2  # noqa: F401

        have_gillespy = True
    except ImportError:
        have_gillespy = False
        print("gillespy2 not installed; using NumPy tau-leap fallback.")

    rows: list[dict] = []
    seed = 0
    cells = [(sw, pr) for sw in switching_rates for pr in pause_rates]
    iterator = track(cells, description="grid cells") if (HAVE_RICH and show_progress) else cells
    for sw, pr in iterator:
        if have_gillespy:
            model = build_model(switching_rate=sw, pause_target_rate=pr)
            model.timespan(range(t_end + 1))
            seed += 1
            # NumPySSASolver: pure Python, runs without a g++ toolchain
            # (SSACSolver is the faster drop-in where g++ exists).
            res = model.run(
                solver=gillespy2.NumPySSASolver,
                number_of_trajectories=n_trajectories,
                seed=seed,
            )
            # run() returns a Results object (len + [i] indexing) on
            # gillespy2>=1.8; older versions returned a plain list.
            # Indexing Results itself (res["A"]) concatenates trajectories --
            # never do that here, it silently collapses n rows into one.
            for i in range(len(res)):
                traj = res[i]
                rows.append(
                    {
                        "switching_rate": sw,
                        "pause_target_rate": pr,
                        "final_active": int(traj["A"][-1]),
                    }
                )
        else:
            for _ in range(n_trajectories):
                seed += 1
                rows.append(
                    {
                        "switching_rate": sw,
                        "pause_target_rate": pr,
                        "final_active": fallback_tauleap(sw, pr, seed=seed),
                    }
                )
    return pd.DataFrame(rows)


def summarize(df: pd.DataFrame) -> pd.DataFrame:
    """Mean final active cells per grid cell (pure; tested)."""
    return (
        df.groupby(["switching_rate", "pause_target_rate"], as_index=False)["final_active"]
        .mean()
        .rename(columns={"final_active": "mean_final_active"})
    )


def render_summary(df: pd.DataFrame) -> None:
    """Pretty-print the grid summary; plain print when rich is missing."""
    summ = summarize(df)
    if not HAVE_RICH:
        print(summ.to_string(index=False))
        return
    table = Table(title="Dormancy grid: mean final active cells")
    table.add_column("switching_rate", justify="right")
    table.add_column("pause_target_rate", justify="right")
    table.add_column("mean_final_active", justify="right")
    for _, r in summ.iterrows():
        table.add_row(
            f"{r['switching_rate']}",
            f"{r['pause_target_rate']}",
            f"{r['mean_final_active']:.1f}",
        )
    Console().print(table)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--switching", nargs="+", type=float, default=[0.005, 0.01, 0.02])
    ap.add_argument("--pause", nargs="+", type=float, default=[0.0, 0.05, 0.1])
    ap.add_argument("--n-trajectories", type=int, default=50)
    ap.add_argument("--t-end", type=int, default=500)
    ap.add_argument("--quiet", action="store_true", help="disable progress bar")
    args = ap.parse_args()

    df = run_grid(
        args.switching,
        args.pause,
        args.n_trajectories,
        args.t_end,
        show_progress=not args.quiet,
    )
    out = ROOT / "results" / "m4"
    out.mkdir(parents=True, exist_ok=True)
    df.to_csv(out / "dormancy_grid.csv", index=False)
    render_summary(df)


if __name__ == "__main__":
    main()
