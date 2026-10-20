
from pathlib import Path
import re
import requests

import pandas as pd
from string_api import map_gene_ids

REFERENCE_FILE = Path(
    "data/reference/reactome_human_protein_references.csv"
)
OUTPUT_DIR = Path("data/reference")


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    references = pd.read_csv(REFERENCE_FILE)

    # Extract unique human candidate symbols.
    genes = sorted(
        references["gene_symbol"].dropna().astype(str).unique()
    )

    # Basic format audit; this is not proof of gene validity.
    symbol_pattern = re.compile(r"^[A-Z][A-Z0-9-]*$")
    flagged = [
        gene for gene in genes
        if len(gene) < 2 or not symbol_pattern.fullmatch(gene)
    ]

    print(f"Human candidate symbols: {len(genes)}")
    print(f"Symbols flagged for review: {len(flagged)}")
    print("Flagged symbols:", ", ".join(flagged) if flagged else "None")

    # Preserve a snapshot of the candidate list.
    candidate_file = (
        OUTPUT_DIR / "interferon_genes_reactome_human_candidates.txt"
    )
    candidate_file.write_text(
        "\n".join(genes) + "\n",
        encoding="utf-8",
    )

    # Map symbols using the existing STRING API function.
    mapping = map_gene_ids(genes)

    if mapping.empty:
        raise ValueError(
            "STRING returned no mappings. Check the API response "
            "before proceeding."
        )

    mapping.to_csv(
        OUTPUT_DIR / "reactome_string_mapping_raw.csv",
        index=False,
    )

    print("\nSTRING mapping columns:")
    print(", ".join(mapping.columns))

    # Identify the query and STRING identifier columns.
    query_col = next(
        (c for c in ["queryItem", "query", "input"] if c in mapping.columns),
        None,
    )
    id_col = next(
        (c for c in ["stringId", "string_id"] if c in mapping.columns),
        None,
    )

    if query_col is None or id_col is None:
        raise ValueError(
            "Could not identify the query and STRING ID columns. "
            "Inspect reactome_string_mapping_raw.csv."
        )

    mapping[query_col] = mapping[query_col].astype(str).str.strip()
    mapping[id_col] = mapping[id_col].astype(str).str.strip()

    # Summarize mapping results for every candidate symbol.
    mapped_ids = (
        mapping.dropna(subset=[id_col])
        .groupby(query_col)[id_col]
        .agg(lambda values: sorted(set(values)))
        .to_dict()
    )

    rows = []

    for gene in genes:
        ids = mapped_ids.get(gene, [])
        if not isinstance(ids, list):
            ids = []

        rows.append({
            "gene_symbol": gene,
            "string_id_count": len(ids),
            "string_ids": ";".join(ids),
            "mapping_status": (
                "unmapped" if len(ids) == 0
                else "mapped_unique" if len(ids) == 1
                else "multiple_matches"
            ),
            "format_flag": gene in flagged,
        })

    audit = pd.DataFrame(rows)
    audit.to_csv(
        OUTPUT_DIR / "reactome_gene_curation_report.csv",
        index=False,
    )
    
    # Export STRING-mappable and unmapped symbols separately.
    mapped_genes = sorted(
        audit.loc[
            audit["mapping_status"] == "mapped_unique",
            "gene_symbol",
        ].unique()
    )

    unmapped_genes = sorted(
        audit.loc[
            audit["mapping_status"] == "unmapped",
            "gene_symbol",
        ].unique()
    )

    mapped_file = (
        OUTPUT_DIR / "interferon_genes_reactome_string_mappable.txt"
    )
    mapped_file.write_text(
        "\n".join(mapped_genes) + "\n",
        encoding="utf-8",
    )

    unmapped_file = (
        OUTPUT_DIR / "interferon_genes_reactome_unmapped.txt"
    )
    unmapped_file.write_text(
        "\n".join(unmapped_genes) + "\n"
        if unmapped_genes else "",
        encoding="utf-8",
    )

    print(f"\nHuman candidates: {len(genes)}")
    print(f"Uniquely STRING-mapped: {len(mapped_genes)}")
    print(f"Unmapped: {len(unmapped_genes)}")
    print(f"Saved: {mapped_file}")
    print(f"Saved: {unmapped_file}")

    print("\nMapping summary:")
    print(audit["mapping_status"].value_counts().to_string())

    unmapped = audit[audit["mapping_status"] == "unmapped"]

    print("\nUnmapped gene symbols:")
    if unmapped.empty:
        print("None")
    else:
        print(unmapped["gene_symbol"].to_string(index=False))
        
    if not unmapped.empty:
        print("\nReactome source records for unmapped symbols:")

        source_records = references[
            references["gene_symbol"].isin(
                unmapped["gene_symbol"]
            )
        ]

        columns = [
            "gene_symbol",
            "reference_display_name",
            "reference_stable_id",
            "reference_schema_class",
            "participant_name",
        ]

        print(
            source_records[columns]
            .drop_duplicates()
            .to_string(index=False)
        )
        
        print("\nUniProt validation for unmapped symbols:")

        for _, row in source_records[
            ["gene_symbol", "reference_display_name"]
        ].drop_duplicates().iterrows():
            match = re.search(
                r"UniProt:([A-Z0-9]+(?:-[0-9]+)?)",
                row["reference_display_name"],
            )

            if not match:
                continue

            accession = match.group(1)

            try:
                response = requests.get(
                    f"https://rest.uniprot.org/uniprotkb/{accession}.json",
                    timeout=30,
                )
                response.raise_for_status()
                data = response.json()

                organism = data.get("organism") or {}
                protein = (
                    data.get("proteinDescription", {})
                    .get("recommendedName", {})
                    .get("fullName", {})
                    .get("value")
                )

                print(
                    f"{row['gene_symbol']} | {accession} | "
                    f"{protein} | {organism.get('scientificName')} | "
                    f"taxon: {organism.get('taxonId')}"
                )

            except requests.RequestException as error:
                print(f"{accession}: validation failed ({error})")

    print("\nSaved:")
    print("data/reference/interferon_genes_reactome_human_candidates.txt")
    print("data/reference/reactome_string_mapping_raw.csv")
    print("data/reference/reactome_gene_curation_report.csv")


if __name__ == "__main__":
    main()
