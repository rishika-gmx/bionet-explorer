
from pathlib import Path
import time

import networkx as nx
import pandas as pd

from string_api import fetch_interactions, read_gene_list, map_gene_ids

RESULTS_DIR = Path("results")
THRESHOLDS = [400, 700, 900]


def build_graph(interactions):
    graph = nx.Graph()

    if interactions.empty:
        return graph

    for _, row in interactions.iterrows():
        a = row["preferredName_A"]
        b = row["preferredName_B"]
        score = float(row["score"])

        graph.add_edge(a, b, weight=score)

    return graph


def analyze_graph(graph, threshold):
    if graph.number_of_nodes() == 0:
        return [], {
            "threshold": threshold / 1000,
            "nodes": 0,
            "edges": 0,
            "density": 0.0,
            "components": 0,
        }

    degree = dict(graph.degree())
    betweenness = nx.betweenness_centrality(graph, normalized=True)

    rankings = pd.DataFrame({
        "protein": list(graph.nodes()),
        "degree": [degree[n] for n in graph.nodes()],
        "betweenness_centrality": [
            betweenness[n] for n in graph.nodes()
        ],
    })

    rankings = rankings.sort_values(
        ["betweenness_centrality", "degree"],
        ascending=False,
    )
    rankings["threshold"] = threshold / 1000

    summary = {
        "threshold": threshold / 1000,
        "nodes": graph.number_of_nodes(),
        "edges": graph.number_of_edges(),
        "density": nx.density(graph),
        "components": nx.number_connected_components(graph),
    }

    return rankings.to_dict("records"), summary


def main():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    genes = read_gene_list()
    mapping = map_gene_ids(genes)
    string_ids = mapping["stringId"].dropna().astype(str).unique().tolist()

    all_rankings = []
    all_summaries = []

    for index, threshold in enumerate(THRESHOLDS):
        print(f"\nFetching interactions at threshold {threshold}...")

        interactions = fetch_interactions(
            string_ids,
            required_score=threshold,
        )

        interactions.to_csv(
            RESULTS_DIR / f"interactions_{threshold}.csv",
            index=False,
        )

        graph = build_graph(interactions)
        rankings, summary = analyze_graph(graph, threshold)

        all_rankings.extend(rankings)
        all_summaries.append(summary)

        print(summary)

        if index < len(THRESHOLDS) - 1:
            time.sleep(1)

    pd.DataFrame(all_summaries).to_csv(
        RESULTS_DIR / "threshold_comparison.csv",
        index=False,
    )

    pd.DataFrame(all_rankings).to_csv(
        RESULTS_DIR / "threshold_centrality.csv",
        index=False,
    )

    print("\nSaved:")
    print("  results/threshold_comparison.csv")
    print("  results/threshold_centrality.csv")
    print("  results/interactions_400.csv")
    print("  results/interactions_700.csv")
    print("  results/interactions_900.csv")


if __name__ == "__main__":
    main()
