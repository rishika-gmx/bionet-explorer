
from io import StringIO
from pathlib import Path
import time

import pandas as pd
import requests

BASE_URL = "https://string-db.org/api"
SPECIES = 9606  # Homo sapiens
CALLER_IDENTITY = "BioNetExplorer"
DATA_DIR = Path("data")
RESULTS_DIR = Path("results")


def read_gene_list(path=DATA_DIR / "interferon_genes.txt"):
    """Read gene symbols, ignoring blank lines, comments and duplicates."""
    with open(path, encoding="utf-8") as file:
        genes = [
            line.strip()
            for line in file
            if line.strip() and not line.lstrip().startswith("#")
        ]

    genes = list(dict.fromkeys(genes))

    if not genes:
        raise ValueError("The input gene list is empty.")

    return genes


def map_gene_ids(genes):
    """Map human gene symbols to STRING identifiers."""
    response = requests.post(
        f"{BASE_URL}/tsv/get_string_ids",
        data={
            "identifiers": "\r".join(genes),
            "species": SPECIES,
            "echo_query": 1,
            "caller_identity": CALLER_IDENTITY,
        },
        timeout=30,
    )
    response.raise_for_status()

    if not response.text.strip():
        raise ValueError("STRING returned no identifier mappings.")

    mapped = pd.read_csv(StringIO(response.text), sep="\t")

    required = {"queryItem", "stringId", "preferredName"}
    missing = required - set(mapped.columns)

    if missing:
        raise ValueError(f"Missing mapping columns: {sorted(missing)}")

    return mapped


def fetch_interactions(string_ids, required_score=700):
    """Retrieve STRING functional associations above the score threshold."""
    if not string_ids:
        raise ValueError("No STRING identifiers were supplied.")

    if not 0 <= required_score <= 1000:
        raise ValueError("Confidence threshold must be between 0 and 1000.")

    response = requests.post(
        f"{BASE_URL}/tsv/network",
        data={
            "identifiers": "\r".join(dict.fromkeys(string_ids)),
            "species": SPECIES,
            "required_score": required_score,
            "network_type": "functional",
            "add_nodes": 0,
            "caller_identity": CALLER_IDENTITY,
        },
        timeout=60,
    )
    response.raise_for_status()

    if not response.text.strip():
        return pd.DataFrame()

    return pd.read_csv(StringIO(response.text), sep="\t")


def main():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    genes = read_gene_list()
    print(f"Input genes: {len(genes)}")

    mapped = map_gene_ids(genes)
    mapped.to_csv(RESULTS_DIR / "gene_mapping.csv", index=False)

    mapped_symbols = set(mapped["queryItem"].astype(str))
    unmapped = [gene for gene in genes if gene not in mapped_symbols]

    print(f"Mapped input symbols: {len(mapped_symbols)}")
    print(f"Unmapped input symbols: {len(unmapped)}")

    if unmapped:
        print("Unmapped genes:", ", ".join(unmapped))

    string_ids = mapped["stringId"].dropna().astype(str).unique().tolist()

    if not string_ids:
        raise ValueError("No valid STRING identifiers were mapped.")

    # Be considerate to the STRING API between consecutive requests.
    time.sleep(1)

    interactions = fetch_interactions(string_ids, required_score=700)
    interactions.to_csv(RESULTS_DIR / "interactions.csv", index=False)

    print(f"Mapped identifiers: {len(string_ids)}")
    print(f"Interaction rows: {len(interactions)}")
    print("Saved results/gene_mapping.csv")
    print("Saved results/interactions.csv")


if __name__ == "__main__":
    main()
    