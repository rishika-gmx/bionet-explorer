
from pathlib import Path

import networkx as nx
import pandas as pd

RESULTS_DIR = Path("results")


def build_network(path=RESULTS_DIR / "interactions.csv"):
    df = pd.read_csv(path)

    required = {
        "stringId_A",
        "stringId_B",
        "preferredName_A",
        "preferredName_B",
        "score",
    }

    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    graph = nx.Graph()

    for _, row in df.iterrows():
        protein_a = row["preferredName_A"]
        protein_b = row["preferredName_B"]
        score = float(row["score"])

        graph.add_node(protein_a, string_id=row["stringId_A"])
        graph.add_node(protein_b, string_id=row["stringId_B"])

        graph.add_edge(
            protein_a,
            protein_b,
            weight=score,
            distance=1.0 / score,
        )

    return graph


def calculate_centrality(graph):
    degree = dict(graph.degree())
    weighted_degree = dict(graph.degree(weight="weight"))

    betweenness_unweighted = nx.betweenness_centrality(
        graph,
        normalized=True,
        weight=None,
    )

    betweenness_weighted = nx.betweenness_centrality(
        graph,
        normalized=True,
        weight="distance",
    )

    results = pd.DataFrame({
        "protein": list(graph.nodes()),
        "degree": [degree[node] for node in graph.nodes()],
        "weighted_degree": [
            weighted_degree[node] for node in graph.nodes()
        ],
        "betweenness_unweighted": [
            betweenness_unweighted[node] for node in graph.nodes()
        ],
        "betweenness_weighted": [
            betweenness_weighted[node] for node in graph.nodes()
        ],
    })

    return results.sort_values(
        by=["betweenness_weighted", "degree"],
        ascending=False,
    )


def main(
    input_path=RESULTS_DIR / "interactions.csv",
    output_dir=RESULTS_DIR,
):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    graph = build_network(input_path)

    print("NETWORK SUMMARY")
    print("-" * 32)
    print(f"Nodes: {graph.number_of_nodes()}")
    print(f"Edges: {graph.number_of_edges()}")
    print(f"Connected components: {nx.number_connected_components(graph)}")
    print(f"Density: {nx.density(graph):.3f}")

    centrality = calculate_centrality(graph)
    output_file = output_dir / "centrality_results.csv"
    centrality.to_csv(output_file, index=False)

    print("\nTOP 10 PROTEINS BY WEIGHTED BETWEENNESS CENTRALITY")
    print(
        centrality.head(10).to_string(
            index=False,
            float_format=lambda value: f"{value:.4f}",
        )
    )

    print(f"\nSaved {output_file}")


if __name__ == "__main__":
    main()
