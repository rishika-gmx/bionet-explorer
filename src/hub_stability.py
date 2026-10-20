
from pathlib import Path

import pandas as pd

RESULTS_DIR = Path("results")
INPUT_FILE = RESULTS_DIR / "threshold_centrality.csv"
OUTPUT_FILE = RESULTS_DIR / "hub_stability.csv"

CANDIDATES = ["STAT1", "STAT2", "IRF9", "IFNAR1", "JAK1"]


def main():
    df = pd.read_csv(INPUT_FILE)

    # Rank proteins by betweenness within each threshold.
    df["betweenness_rank"] = (
        df.groupby("threshold")["betweenness_centrality"]
        .rank(method="min", ascending=False)
    )

    # Also calculate degree ranks, preserving ties.
    df["degree_rank"] = (
        df.groupby("threshold")["degree"]
        .rank(method="min", ascending=False)
    )

    subset = df[df["protein"].isin(CANDIDATES)].copy()

    summary = subset.pivot(
        index="protein",
        columns="threshold",
        values="betweenness_rank",
    )

    summary.columns = [
        f"betweenness_rank_{threshold:.1f}"
        for threshold in summary.columns
    ]

    summary.to_csv(OUTPUT_FILE)

    print("Betweenness ranks across thresholds:")
    print(summary.to_string())

    print(f"\nSaved: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
