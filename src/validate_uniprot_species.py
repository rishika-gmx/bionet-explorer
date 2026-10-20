
from pathlib import Path
import time

import pandas as pd
import requests

INPUT_FILE = Path("data/reference/reactome_protein_references.csv")
OUTPUT_DIR = Path("data/reference")
HUMAN_TAXON_ID = 9606
UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/{}.json"



def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    references = pd.read_csv(INPUT_FILE)

    references["uniprot_accession"] = references[
        "reference_display_name"
    ].str.extract(r"UniProt:([A-Z0-9]+(?:-[0-9]+)?)")

    accessions = sorted(
        references["uniprot_accession"].dropna().unique()
    )

    print(f"Unique UniProt accessions to check: {len(accessions)}")

    session = requests.Session()
    session.headers.update({"Accept": "application/json"})

    results = []

    for index, accession in enumerate(accessions, start=1):
        try:
            response = session.get(
                UNIPROT_URL.format(accession),
                timeout=30,
            )

            if response.status_code == 429 or response.status_code >= 500:
                time.sleep(2)
                response = session.get(
                    UNIPROT_URL.format(accession),
                    timeout=30,
                )

            if response.status_code == 404:
                results.append({
                    "uniprot_accession": accession,
                    "organism": None,
                    "taxon_id": None,
                    "status": "not_found",
                })
                continue

            response.raise_for_status()
            data = response.json()
            organism = data.get("organism") or {}
            taxon_id = organism.get("taxonId")

            results.append({
                "uniprot_accession": accession,
                "organism": organism.get("scientificName"),
                "taxon_id": taxon_id,
                "status": (
                    "human"
                    if taxon_id == HUMAN_TAXON_ID
                    else "nonhuman"
                ),
            })

        except requests.RequestException as error:
            results.append({
                "uniprot_accession": accession,
                "organism": None,
                "taxon_id": None,
                "status": f"error: {type(error).__name__}",
            })

        if index % 25 == 0 or index == len(accessions):
            print(f"Checked {index}/{len(accessions)} accessions")

        time.sleep(0.1)

    report = pd.DataFrame(results)
    report_file = OUTPUT_DIR / "uniprot_species_validation.csv"
    report.to_csv(report_file, index=False)

    print("\nSpecies validation summary:")
    print(report["status"].value_counts(dropna=False).to_string())

    # Join validation results to Reactome reference records.
    audited = references.merge(
        report,
        on="uniprot_accession",
        how="left",
    )

    nonhuman = audited[audited["status"] == "nonhuman"].copy()
    nonhuman_file = OUTPUT_DIR / "reactome_nonhuman_references.csv"
    nonhuman.to_csv(nonhuman_file, index=False)

    if nonhuman.empty:
        print("\nNo nonhuman references found.")
    else:
        columns = [
            "gene_symbol",
            "reference_display_name",
            "organism",
            "taxon_id",
            "participant_name",
        ]

        print("\nNonhuman Reactome references:")
        print(
            nonhuman[columns]
            .drop_duplicates()
            .sort_values(["organism", "gene_symbol"])
            .to_string(index=False)
        )

    
    # Export only references validated as human.
    human_references = audited[audited["status"] == "human"].copy()

    human_reference_file = (
        OUTPUT_DIR / "reactome_human_protein_references.csv"
    )
    human_references.to_csv(human_reference_file, index=False)

    # Create a deduplicated candidate gene-symbol list.
    human_genes = sorted(
        human_references["gene_symbol"].dropna().unique()
    )

    human_gene_file = (
        OUTPUT_DIR / "interferon_genes_reactome_human_candidates.txt"
    )
    human_gene_file.write_text(
        "\n".join(human_genes) + "\n",
        encoding="utf-8",
    )

    unresolved = audited[
        audited["status"].astype(str).str.startswith("error:")
        | audited["status"].isin(["not_found"])
        | audited["status"].isna()
    ]

    print("\nHuman candidate summary:")
    print(f"Human protein-reference rows: {len(human_references)}")
    print(f"Unique human candidate symbols: {len(human_genes)}")
    print(f"Unresolved reference rows: {len(unresolved)}")

    print(f"\nSaved: {report_file}")
    print(f"Saved: {nonhuman_file}")
    print(f"Saved: {human_reference_file}")
    print(f"Saved: {human_gene_file}")


if __name__ == "__main__":
    main()
