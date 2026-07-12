#!/usr/bin/env python3
"""Create a small deterministic source-recipient frame for smoke testing."""

from pathlib import Path

import numpy as np
import pandas as pd


def main() -> None:
    rng = np.random.default_rng(20260711)
    rows = []
    for source in range(60):
        treat = int(source >= 30)
        mouse = f"M{1 + source % 3}"
        source_effect = rng.normal(0, 0.45)
        for recipient in range(rng.integers(80, 150)):
            distance = rng.uniform(2, 250)
            near = np.exp(-distance / 35)
            noise = rng.normal(0, 0.55)
            rows.append(
                {
                    "source_id": f"S{source:03d}",
                    "recipient_id": f"S{source:03d}_R{recipient:03d}",
                    "mouse": mouse,
                    "treat": treat,
                    "dist_um": distance,
                    "log_depth": rng.normal(6.8, 0.35),
                    "program_a": 0.16 * treat * near + source_effect + noise,
                    "program_b": source_effect + noise,
                }
            )
    out = Path("data/synthetic/source_recipient_frame.csv")
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out, index=False)
    print(out)


if __name__ == "__main__":
    main()

