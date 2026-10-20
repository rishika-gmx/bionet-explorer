
from pathlib import Path
import pandas as pd

from string_api import map_gene_ids, fetch_interactions

GENE_FILE = Path(
    "data/reference/interferon_genes_reactome_string_mappable.txt"
)
OUTPUT_DIR = Path("results/reference_network")


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    genes = [
        line.strip()
        for line in GENE_FILE.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]

    print(f"Input genes: {len(genes)}")

    mapping = map_gene_ids(genes)
    if mapping.empty:
        raise ValueError("STRING returned no gene mappings.")

    mapping_file = OUTPUT_DIR / "gene_mapping.csv"
    mapping.to_csv(mapping_file, index=False)

    mapped_ids = (
        mapping.dropna(subset=["stringId"])["stringId"]
        .drop_duplicates()
        .tolist()
    )

    print(f"Unique STRING protein IDs: {len(mapped_ids)}")

    interactions = fetch_interactions(
        mapped_ids,
        required_score=700,
    )

    if interactions.empty:
        raise ValueError(
            "STRING returned no interactions. Check the API response."
        )

    interaction_file = OUTPUT_DIR / "interactions.csv"
    interactions.to_csv(interaction_file, index=False)

    print(f"Interaction rows: {len(interactions)}")
    print(f"Saved: {mapping_file}")
    print(f"Saved: {interaction_file}")


if __name__ == "__main__":
    main()
