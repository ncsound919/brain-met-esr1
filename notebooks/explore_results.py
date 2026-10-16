"""Explore brain-met-esr1 results interactively.

Run:  marimo edit notebooks/explore_results.py   (needs `pip install brain-met-esr1[explore]`)
Or:   python notebooks/explore_results.py        (executes cells as a script)

Reactive: moving the lambda slider re-runs only the attenuation cells.
"""

import marimo

app = marimo.App()


@app.cell
def __():
    import marimo as mo
    from pathlib import Path

    ROOT = Path(__file__).resolve().parents[1]
    return (mo, ROOT)


@app.cell
def __(mo):
    lam_slider = mo.ui.slider(0.5, 0.95, step=0.05, value=0.7, label="λ (reliability)")
    lam_slider
    return (lam_slider,)


@app.cell
def __(lam_slider, mo, ROOT):
    import math

    import numpy as np
    import pandas as pd

    lam = lam_slider.value
    # Pooled anchor: recomputed from study_table.parquet (Woolf SE);
    # meta anchor: external 0.35 (0.24-0.52), abstract-only.
    study_path = ROOT / "data" / "raw" / "study_table.parquet"
    triples = [("meta", 0.35, 0.24, 0.52)]
    if study_path.exists():
        st = pd.read_parquet(study_path)
        a = st["esr1_brain_mut"] + 0.5
        b = st["n_brain"] - st["esr1_brain_mut"] + 0.5
        c = st["esr1_nonCNS_mut"] + 0.5
        d = st["n_nonCNS"] - st["esr1_nonCNS_mut"] + 0.5
        lori = np.log(a * d / (b * c))
        sei = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
        w = 1 / sei**2  # fixed-effect inverse-variance pooling (display only)
        lor = float((w * lori).sum() / w.sum())
        se = float(np.sqrt(1 / w.sum()))
        triples.append(
            ("pooled", math.exp(lor), math.exp(lor - 1.96 * se), math.exp(lor + 1.96 * se))
        )
    corrected_df = pd.DataFrame(
        {
            "anchor": [t[0] for t in triples for _ in range(3)],
            "estimate": ["OR_true", "lo", "hi"] * len(triples),
            "value": [math.exp(math.log(v) / lam) for t in triples for v in t[1:]],
        }
    )
    mo.md(f"## M1 attenuation at λ = {lam}\n`OR_true = exp(log(OR_obs) / λ)`")
    return (corrected_df, lam)


@app.cell
def __(corrected_df):
    corrected_df
    return


@app.cell
def __(ROOT, mo):
    import pandas as pd

    study_path = ROOT / "data" / "raw" / "study_table.parquet"
    if study_path.exists():
        study_df = pd.read_parquet(study_path)
        mo.md(f"## Study table ({len(study_df)} studies)")
    else:
        study_df = pd.DataFrame()
        mo.md("## Study table\n`study_table.parquet` missing — run `make_synthetic.py`.")
    return (study_df,)


@app.cell
def __(study_df):
    study_df
    return


@app.cell
def __(ROOT, mo):
    import pandas as pd

    dorm_path = ROOT / "results" / "m4" / "dormancy_grid.csv"
    if dorm_path.exists():
        dorm_df = (
            pd.read_csv(dorm_path)
            .groupby(["switching_rate", "pause_target_rate"], as_index=False)["final_active"]
            .mean()
        )
        mo.md("## M4 dormancy grid means")
    else:
        dorm_df = pd.DataFrame()
        mo.md("## M4 dormancy grid\nNo `results/m4/dormancy_grid.csv` yet — run `simulate.py`.")
    return (dorm_df,)


@app.cell
def __(dorm_df):
    dorm_df
    return


if __name__ == "__main__":
    app.run()
