#!/usr/bin/env python3
"""Run the source-aware audit on a long-form analysis frame."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from perturbation_psf.audit import (
    benjamini_hochberg,
    compare_replication_units,
    equal_source_test,
    wild_cluster_bootstrap,
)


def read_frame(path: Path) -> pd.DataFrame:
    if path.suffix == ".parquet":
        return pd.read_parquet(path)
    return pd.read_csv(path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, default=Path("data/generated"))
    parser.add_argument("--outcomes", nargs="+", required=True)
    parser.add_argument("--covariates", nargs="*", default=[])
    parser.add_argument("--wild-reps", type=int, default=999)
    args = parser.parse_args()

    frame = read_frame(args.input)
    args.outdir.mkdir(parents=True, exist_ok=True)
    rows, wild, equal = [], [], []
    for outcome in args.outcomes:
        rows.append(
            compare_replication_units(
                frame, outcome, covariates=args.covariates
            ).to_dict()
        )
        wild.append(
            wild_cluster_bootstrap(
                frame,
                outcome,
                covariates=args.covariates,
                reps=args.wild_reps,
            )
        )
        equal.append(equal_source_test(frame, outcome, near_field="dist_um <= 40"))

    audit = pd.DataFrame(rows)
    audit["q_naive"] = benjamini_hochberg(audit.p_naive)
    audit["q_source"] = benjamini_hochberg(audit.p_source)
    wild_df = pd.DataFrame(wild)
    wild_df["q_wild"] = benjamini_hochberg(wild_df.p_wild)
    audit.to_csv(args.outdir / "replication_unit_audit.csv", index=False)
    wild_df.to_csv(args.outdir / "wild_cluster_bootstrap.csv", index=False)
    pd.DataFrame(equal).to_csv(args.outdir / "equal_source.csv", index=False)

    summary = {
        "outcomes": len(args.outcomes),
        "naive_q05": int((audit.q_naive <= 0.05).sum()),
        "source_q05": int((audit.q_source <= 0.05).sum()),
        "wild_q05": int((wild_df.q_wild <= 0.05).sum()),
        "median_se_inflation": float(audit.se_inflation.median()),
    }
    (args.outdir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

