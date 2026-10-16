"""Provenance gate: a table is not promotable until something declared it real.

Why this exists
---------------
`data/make_synthetic.py` writes schema-valid placeholders so the pipeline runs
before extraction lands. Three of them encode the hypothesis by construction:

  * `first_progression`  -- `h_brain = 0.001 if esr1=="mut" else 0.003`
  * `expression`         -- `shift = 1.5 if esr1=="mut"` (the *falsifying* case)

So M5's competing-risks output is a restatement of a hard-coded hazard ratio,
not evidence. It once sat in `results/m5/` behind a `print()`-and-continue.

Detecting synthetic data from the data alone is not sufficient
-------------------------------------------------------------
`detection_history` and `expression_samples` carry `assay_platform ==
"synthetic"`, so a column scan catches them. **`first_progression` has no
`assay_platform` column at all** (see schema.py::FirstProgression), and neither
does the tidy `expression` matrix. A scan keyed on any column therefore passes
the two tables that matter most. Provenance is consequently recorded by the
*writer*, in a manifest, and an undeclared table is `unknown` -- which is
treated as not-real. Missing provenance fails closed.

CLI:
    python data/provenance.py status
    python data/provenance.py require first_progression.parquet expression.parquet
    python data/provenance.py stamp --table first_progression.parquet --status real --source "<url/doi>"
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Literal

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
MANIFEST = RAW / "provenance.json"

Status = Literal["real", "synthetic", "unknown"]

#: Statuses that may back a published finding. `unknown` is deliberately absent:
#: an undeclared table is never promotable.
PROMOTABLE: frozenset[str] = frozenset({"real"})

MANIFEST_VERSION = 1


class ProvenanceError(RuntimeError):
    """Raised when a result would rest on data not known to be real."""


def _read() -> dict:
    if not MANIFEST.exists():
        return {"version": MANIFEST_VERSION, "tables": {}}
    try:
        data = json.loads(MANIFEST.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise ProvenanceError(f"provenance manifest is unreadable: {MANIFEST} ({exc})") from exc
    if not isinstance(data.get("tables"), dict):
        raise ProvenanceError(f"provenance manifest has no 'tables' object: {MANIFEST}")
    return data


def stamp(table: str, status: Status, source: str, note: str = "") -> None:
    """Declare where a raw table came from. Called by the writer, not by a reader."""
    if status not in ("real", "synthetic"):
        raise ValueError(f"status must be 'real' or 'synthetic', got {status!r}")
    if not source.strip():
        raise ValueError(f"refusing to stamp {table!r} with an empty source")
    data = _read()
    data["version"] = MANIFEST_VERSION
    data["tables"][table] = {
        "status": status,
        "source": source,
        "note": note,
        "declared_at": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "declared_by": Path(sys.argv[0]).name if sys.argv and sys.argv[0] else "unknown",
    }
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    MANIFEST.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def table_status(table: str) -> Status:
    entry = _read()["tables"].get(table)
    if not entry:
        return "unknown"
    status = entry.get("status")
    return status if status in ("real", "synthetic") else "unknown"


def describe(table: str) -> str:
    entry = _read()["tables"].get(table)
    if not entry:
        return f"{table}: UNKNOWN (never declared; not promotable)"
    return f"{table}: {entry['status']} <- {entry['source']}"


def frame_is_synthetic(df) -> bool:
    """Belt-and-braces: catch a synthetic frame that was never stamped.

    Only works for tables carrying `assay_platform`. It is a *secondary* check
    by design -- see the module docstring for why the manifest is the primary
    one. Returns False for tables with no such column, which is not an all-clear.
    """
    if "assay_platform" not in getattr(df, "columns", []):
        return False
    return bool(df["assay_platform"].astype(str).str.lower().eq("synthetic").any())


def require_real(tables, allow_synthetic: bool = False) -> dict[str, str]:
    """Fail closed unless every named table is declared real.

    Returns {table: source} on success so callers can record provenance in the
    artefact they emit. With `allow_synthetic` the call succeeds but the returned
    statuses mark the run as non-promotable; callers must stamp their output.
    """
    statuses = {t: table_status(t) for t in tables}
    blocking = {t: s for t, s in statuses.items() if s not in PROMOTABLE}
    if blocking and not allow_synthetic:
        detail = "; ".join(
            f"{t}={s}" + (" (never declared)" if s == "unknown" else "") for t, s in blocking.items()
        )
        raise ProvenanceError(
            f"refusing to emit a result on non-real data: {detail}. "
            "These modules currently run on synthetic placeholders that encode the "
            "hypothesis by construction. Pass --allow-synthetic for a plumbing check "
            "only; such output is never promotable as a finding."
        )
    return statuses


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)

    sub.add_parser("status", help="print every declared table")

    req = sub.add_parser("require", help="exit non-zero unless all tables are real")
    req.add_argument("tables", nargs="+")
    req.add_argument("--allow-synthetic", action="store_true")

    st = sub.add_parser("stamp", help="declare a table's provenance")
    st.add_argument("--table", required=True)
    st.add_argument("--status", required=True, choices=["real", "synthetic"])
    st.add_argument("--source", required=True)
    st.add_argument("--note", default="")

    args = ap.parse_args(argv)

    if args.cmd == "status":
        data = _read()["tables"]
        if not data:
            print("provenance manifest absent: every table is UNKNOWN (not promotable)")
            return 1
        for table in sorted(data):
            print(describe(table))
        return 0

    if args.cmd == "require":
        try:
            statuses = require_real(args.tables, allow_synthetic=args.allow_synthetic)
        except ProvenanceError as exc:
            print(f"PROVENANCE GATE: {exc}", file=sys.stderr)
            return 2
        for table, status in statuses.items():
            print(f"OK {table}: {status}")
        return 0

    stamp(args.table, args.status, args.source, args.note)
    print(f"stamped {args.table} = {args.status}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())