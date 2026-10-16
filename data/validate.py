"""Validate all Parquet tables against data/schema.py. Used by Snakefile."""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from data.schema import (
    DetectionHistory,
    ExpressionRow,
    ExpressionSample,
    FirstProgression,
    StudyRecord,
    check_duplicates,
)

REQUIRED = [
    "detection_history.parquet",
    "study_table.parquet",
    "first_progression.parquet",
    "expression_samples.parquet",
    "expression.parquet",
]


def main() -> None:
    paths = {name: ROOT / "data" / "raw" / name for name in REQUIRED}
    for name, path in paths.items():
        if not path.exists():
            hint = (
                "`python data/ingest_ecosystem.py`"
                if name in ("study_table.parquet", "meta_anchor.json")
                else "`python data/make_synthetic.py`"
            )
            raise SystemExit(f"Missing {path}. Run {hint}.")
    tables = {name: pd.read_parquet(path) for name, path in paths.items()}

    # Parquet round-trips None -> NaN in float columns, and float-dtype
    # .where(..., None) keeps NaN. Sanitize per-record instead.
    def clean(rows: list[dict]) -> list[dict]:
        out = []
        for r in rows:
            out.append(
                {k: (None if isinstance(v, float) and pd.isna(v) else v) for k, v in r.items()}
            )
        return out

    det_rows = clean(tables["detection_history.parquet"].to_dict("records"))
    study_rows = clean(tables["study_table.parquet"].to_dict("records"))
    prog_rows = clean(tables["first_progression.parquet"].to_dict("records"))
    es_rows = clean(tables["expression_samples.parquet"].to_dict("records"))
    ex_rows = clean(tables["expression.parquet"].to_dict("records"))

    check_duplicates(det_rows)
    for row in det_rows:
        DetectionHistory(**row)
    for row in study_rows:
        StudyRecord(**row)
    seen_patients: set[str] = set()
    for row in prog_rows:
        rec = FirstProgression(**row)
        if rec.patient_id in seen_patients:
            raise ValueError(f"Duplicate first-progression row for {rec.patient_id}")
        seen_patients.add(rec.patient_id)
    sample_ids: set[str] = set()
    for row in es_rows:
        rec = ExpressionSample(**row)
        if rec.sample_id in sample_ids:
            raise ValueError(f"Duplicate expression sample {rec.sample_id}")
        sample_ids.add(rec.sample_id)
    for row in ex_rows:
        rec = ExpressionRow(**row)
        if rec.sample_id not in sample_ids:
            raise ValueError(f"Expression row for unknown sample {rec.sample_id}")
    ex_keys = [(r["sample_id"], r["gene"]) for r in ex_rows]
    if len(set(ex_keys)) != len(ex_keys):
        raise ValueError("Duplicate (sample_id, gene) expression rows: pivot would fail")

    stamp = ROOT / "data" / "validated" / ".stamp"
    stamp.parent.mkdir(parents=True, exist_ok=True)
    stamp.write_text("ok\n")
    print(
        f"Validated {len(det_rows)} detection, {len(study_rows)} study, "
        f"{len(prog_rows)} progression, {len(es_rows)} samples, "
        f"{len(ex_rows)} expression rows."
    )


if __name__ == "__main__":
    main()
