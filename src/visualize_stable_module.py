
from pathlib import Path

import networkx as nx
import pandas as pd
import plotly.graph_objects as go

from network_analysis import build_network

RESULTS_DIR = Path("results/reference_network")
INTERACTION_FILE = RESULTS_DIR / "interactions.csv"
COMMUNITY_FILE = RESULTS_DIR / "community_assignments.csv"
CENTRALITY_FILE = RESULTS_DIR / "centrality_results.csv"
OUTPUT_FILE = RESULTS_DIR / "stable_interferon_module.html"


def main():
    assignments = pd.read_csv(COMMUNITY_FILE)

    low = assignments[assignments["threshold"] == 0.7]
    high = assignments[assignments["threshold"] == 0.9]

    low_community = set(low.loc[low["community"] == 4, "protein"])
    high_community = set(high.loc[high["community"] == 8, "protein"])

    stable_proteins = low_community & high_community

    graph = build_network(INTERACTION_FILE)
    module = graph.subgraph(stable_proteins).copy()

    centrality = pd.read_csv(CENTRALITY_FILE).set_index("protein")

    positions = nx.spring_layout(
        module,
        seed=42,
        weight="weight",
        k=0.8,
        iterations=150,
    )

    edge_x, edge_y = [], []

    for source, target in module.edges():
        x0, y0 = positions[source]
        x1, y1 = positions[target]
        edge_x.extend([x0, x1, None])
        edge_y.extend([y0, y1, None])

    edge_trace = go.Scatter(
        x=edge_x,
        y=edge_y,
        mode="lines",
        line=dict(width=0.7, color="rgba(120,120,120,0.4)"),
        hoverinfo="skip",
        showlegend=False,
    )

    nodes = sorted(module.nodes())
    node_x = [positions[node][0] for node in nodes]
    node_y = [positions[node][1] for node in nodes]

    degrees = [module.degree(node) for node in nodes]
    betweenness = [
        centrality.loc[node, "betweenness_weighted"]
        for node in nodes
    ]

    hover_text = [
        (
            f"<b>{node}</b>"
            f"<br>Module degree: {module.degree(node)}"
            f"<br>Full-network degree: {centrality.loc[node, 'degree']}"
            f"<br>Weighted betweenness: {centrality.loc[node, 'betweenness_weighted']:.4f}"
        )
        for node in nodes
    ]

    node_trace = go.Scatter(
        x=node_x,
        y=node_y,
        mode="markers+text",
        text=nodes,
        textposition="top center",
        textfont=dict(size=9),
        hovertext=hover_text,
        hoverinfo="text",
        marker=dict(
            size=[8 + min(degree, 20) * 0.6 for degree in degrees],
            color=betweenness,
            colorscale="Viridis",
            showscale=True,
            colorbar=dict(title="Full-network<br>weighted betweenness"),
            line=dict(width=0.7, color="white"),
        ),
        showlegend=False,
    )

    figure = go.Figure(data=[edge_trace, node_trace])
    figure.update_layout(
        title=(
            "BioNet Explorer: Stable Interferon/Antiviral Module"
            "<br><sup>Shared community members at STRING thresholds "
            "0.7 and 0.9 | Node size: module degree</sup>"
        ),
        template="plotly_white",
        hovermode="closest",
        showlegend=False,
        margin=dict(l=20, r=20, t=90, b=20),
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
    )

    figure.write_html(str(OUTPUT_FILE), include_plotlyjs=True)

    print(f"Saved: {OUTPUT_FILE}")
    print(f"Shared proteins selected: {len(stable_proteins)}")
    print(f"Edges among selected proteins: {module.number_of_edges()}")
    print(f"Connected components in focused view: {nx.number_connected_components(module)}")


if __name__ == "__main__":
    main()
