
from pathlib import Path
import re
import requests
import pandas as pd

PATHWAY_ID = "R-HSA-913531"

API_URL = (
    "https://reactome.org/ContentService/data/"
    f"participants/{PATHWAY_ID}"
)

OUTPUT_DIR = Path("data/reference")


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    response = requests.get(API_URL, timeout=60)
    response.raise_for_status()
    records = response.json()

    raw_file = OUTPUT_DIR / "reactome_participants.json"
    raw_file.write_text(response.text, encoding="utf-8")

    rows = []

    for participant in records:
        for reference in participant.get("refEntities") or []:
            # Keep protein reference entities, not small molecules,
            # DNA/RNA sequences, or drug references.
            if reference.get("schemaClass") not in {
                "ReferenceGeneProduct",
                "ReferenceIsoform",
            }:
                continue

            display_name = reference.get("displayName", "")

            # Reactome's UniProt reference display names commonly end
            # with a gene symbol, e.g. "UniProt:P42224 STAT1".
            match = re.search(
                r"UniProt:[^\s]+\s+([A-Z0-9-]+)$",
                display_name,
            )

            if not match:
                continue

            rows.append({
                "gene_symbol": match.group(1),
                "reference_display_name": display_name,
                "reference_stable_id": reference.get("stId"),
                "reference_schema_class": reference.get("schemaClass"),
                "participant_name": participant.get("displayName"),
                "participant_class": participant.get("schemaClass"),
            })

    entities = pd.DataFrame(rows)

    if entities.empty:
        raise ValueError(
            "No protein gene symbols extracted. "
            "Inspect the saved Reactome participant JSON."
        )

    entities = entities.drop_duplicates()

    entities.to_csv(
        OUTPUT_DIR / "reactome_protein_references.csv",
        index=False,
    )

    genes = sorted(entities["gene_symbol"].unique())

    gene_file = OUTPUT_DIR / "interferon_genes_reactome_candidates.txt"
    gene_file.write_text("\n".join(genes) + "\n", encoding="utf-8")

    print(f"Pathway participant records: {len(records)}")
    print(f"Protein reference rows: {len(entities)}")
    print(f"Unique candidate gene symbols: {len(genes)}")
    print(f"\nSaved: {raw_file}")
    print(f"Saved: {OUTPUT_DIR / 'reactome_protein_references.csv'}")
    print(f"Saved: {gene_file}")
    print("\nFirst 50 candidate symbols:")
    print(", ".join(genes[:50]))


if __name__ == "__main__":
    main()
