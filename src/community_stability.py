
from itertools import combinations
from pathlib import Path
from math import comb

import pandas as pd

RESULTS_DIR = Path("results/reference_network")
INPUT_FILE = RESULTS_DIR / "community_assignments.csv"
OUTPUT_FILE = RESULTS_DIR / "community_stability.csv"


def adjusted_rand_index(labels_a, labels_b):
    table = pd.crosstab(labels_a, labels_b)

    index = sum(comb(int(n), 2) for n in table.to_numpy().ravel())
    row_pairs = sum(comb(int(n), 2) for n in table.sum(axis=1))
    col_pairs = sum(comb(int(n), 2) for n in table.sum(axis=0))
    total_pairs = comb(len(labels_a), 2)

    expected = row_pairs * col_pairs / total_pairs
    maximum = (row_pairs + col_pairs) / 2

    if maximum == expected:
        return 1.0 if index == maximum else 0.0

    return (index - expected) / (maximum - expected)


def main():
    df = pd.read_csv(INPUT_FILE)

    assignments = df.pivot(
        index="protein",
        columns="threshold",
        values="community",
    )

    thresholds = sorted(assignments.columns)

    if len(thresholds) < 2:
        raise ValueError("At least two thresholds are needed for comparison.")

    rows = []

    for threshold_a, threshold_b in combinations(thresholds, 2):
        shared = assignments[[threshold_a, threshold_b]].dropna()

        ari = adjusted_rand_index(
            shared[threshold_a].tolist(),
            shared[threshold_b].tolist(),
        )

        rows.append({
            "threshold_a": threshold_a,
            "threshold_b": threshold_b,
            "shared_proteins": len(shared),
            "adjusted_rand_index": ari,
        })

    report = pd.DataFrame(rows)
    report.to_csv(OUTPUT_FILE, index=False)

    print("COMMUNITY PARTITION STABILITY")
    print("-" * 36)
    print(report.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
    print("\nARI interpretation: 1.0 means identical partitions; values near 0")
    print("indicate agreement close to chance; negative values indicate less.")
    print(f"Saved: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()
