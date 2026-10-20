
from pathlib import Path

import pandas as pd

RESULTS_DIR = Path("results")


def validate_interactions(
    path=RESULTS_DIR / "interactions.csv",
    report_path=RESULTS_DIR / "validation_report.csv",
):
    df = pd.read_csv(path)

    required = {
        "stringId_A",
        "stringId_B",
        "preferredName_A",
        "preferredName_B",
        "ncbiTaxonId",
        "score",
    }

    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    if df.empty:
        raise ValueError("The interaction table is empty.")

    # Check for missing endpoints and invalid confidence scores.
    endpoint_cols = [
        "stringId_A",
        "stringId_B",
        "preferredName_A",
        "preferredName_B",
    ]

    missing_endpoints = int(df[endpoint_cols].isna().any(axis=1).sum())
    scores = pd.to_numeric(df["score"], errors="coerce")

    invalid_scores = int(
        (scores.isna() | (scores < 0) | (scores > 1)).sum()
    )

    nonhuman_rows = int((df["ncbiTaxonId"] != 9606).sum())

    # STRING reports each association as an endpoint pair.
    edge_pairs = df.apply(
        lambda row: tuple(sorted((
            str(row["stringId_A"]),
            str(row["stringId_B"]),
        ))),
        axis=1,
    )

    duplicate_edges = int(edge_pairs.duplicated().sum())
    self_loops = int(
        (df["stringId_A"] == df["stringId_B"]).sum()
    )

    report = {
        "interaction_rows": len(df),
        "unique_proteins": len(
            set(df["stringId_A"]) | set(df["stringId_B"])
        ),
        "missing_endpoint_rows": missing_endpoints,
        "invalid_score_rows": invalid_scores,
        "nonhuman_rows": nonhuman_rows,
        "duplicate_edges": duplicate_edges,
        "self_loops": self_loops,
        "minimum_score": scores.min(),
        "maximum_score": scores.max(),
    }

    print("\nINTERACTION VALIDATION REPORT")
    print("-" * 32)

    for key, value in report.items():
        print(f"{key}: {value}")

        report_path.parent.mkdir(parents=True, exist_ok=True)

    pd.DataFrame([report]).to_csv(
        report_path,
        index=False,
    )

    print(f"\nSaved {report_path}")


if __name__ == "__main__":
    validate_interactions()
