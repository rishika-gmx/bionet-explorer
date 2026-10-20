
from pathlib import Path

import networkx as nx
import pandas as pd
from networkx.algorithms.community import louvain_communities
from networkx.algorithms.community import modularity

RESULTS_DIR = Path("results/reference_network")
THRESHOLDS = [700, 900]


def build_graph(path, threshold):
    df = pd.read_csv(path)
    df = df[pd.to_numeric(df["score"], errors="coerce") >= threshold / 1000]

    graph = nx.Graph()

    for _, row in df.iterrows():
        graph.add_edge(
            row["preferredName_A"],
            row["preferredName_B"],
            weight=float(row["score"]),
        )

    return graph


def analyze_communities(graph, threshold):
    if graph.number_of_nodes() == 0:
        print(f"No nodes found at threshold {threshold / 1000:.1f}.")
        return []

    communities = louvain_communities(
        graph,
        weight="weight",
        resolution=1,
        seed=42,
    )

    score = modularity(graph, communities, weight="weight")
    sorted_communities = sorted(communities, key=lambda group: min(group))

    rows = []

    for community_id, members in enumerate(sorted_communities, start=1):
        for protein in sorted(members):
            rows.append({
                "threshold": threshold / 1000,
                "community": community_id,
                "protein": protein,
                "community_size": len(members),
                "modularity": score,
            })

    print(f"\nThreshold: {threshold / 1000:.1f}")
    print(f"Nodes: {graph.number_of_nodes()}")
    print(f"Edges: {graph.number_of_edges()}")
    print(f"Communities detected: {len(communities)}")
    print(f"Modularity: {score:.4f}")

    for community_id, members in enumerate(sorted_communities, start=1):
        print(
            f"  Community {community_id} ({len(members)} proteins): "
            f"{', '.join(sorted(members))}"
        )

    return rows


def main():
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    path = RESULTS_DIR / "interactions.csv"

    if not path.exists():
        raise FileNotFoundError(f"Missing curated interaction file: {path}")

    all_rows = []

    for threshold in THRESHOLDS:
        graph = build_graph(path, threshold)
        all_rows.extend(analyze_communities(graph, threshold))

    output = RESULTS_DIR / "community_assignments.csv"
    pd.DataFrame(all_rows).to_csv(output, index=False)

    print(f"\nSaved: {output}")


if __name__ == "__main__":
    main()
