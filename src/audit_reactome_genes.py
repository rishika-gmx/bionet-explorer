
from pathlib import Path
import re
import pandas as pd

INPUT_FILE = Path(
    "data/reference/reactome_protein_references.csv"
)
OUTPUT_DIR = Path("data/reference")


def main():
    df = pd.read_csv(INPUT_FILE)

    genes = sorted(df["gene_symbol"].dropna().unique())

    # Flag symbols that need manual review. This is a review flag,
    # not a claim that the symbol is invalid.
    pattern = re.compile(r"^[A-Z][A-Z0-9-]*$")

    flagged = [
        gene for gene in genes
        if not pattern.fullmatch(str(gene))
        or len(str(gene)) < 2
    ]

    print(f"Unique candidate symbols: {len(genes)}")
    print(f"Symbols flagged for review: {len(flagged)}")
    print("\nFlagged symbols:")
    print(", ".join(flagged) if flagged else "None")

    # Count how often each gene appears across participant references.
    counts = (
        df.groupby("gene_symbol")
        .agg(
            reference_rows=("gene_symbol", "size"),
            participant_count=("participant_name", "nunique"),
        )
        .sort_values("participant_count", ascending=False)
    )

    counts.to_csv(OUTPUT_DIR / "reactome_gene_reference_counts.csv")

    print("\nTop 15 genes by number of distinct participant entities:")
    print(counts.head(15).to_string())

    print(
        "\nSaved: "
        "data/reference/reactome_gene_reference_counts.csv"
    )


if __name__ == "__main__":
    main()
