
from pathlib import Path

import networkx as nx
import pandas as pd
import plotly.graph_objects as go

from network_analysis import build_network

RESULTS_DIR = Path("results/reference_network")
INTERACTION_FILE = RESULTS_DIR / "interactions.csv"
CENTRALITY_FILE = RESULTS_DIR / "centrality_results.csv"
COMMUNITY_FILE = RESULTS_DIR / "community_assignments.csv"
OUTPUT_FILE = RESULTS_DIR / "interferon_network.html"

THRESHOLD = 0.7


def main():
    graph = build_network(INTERACTION_FILE)

    centrality = pd.read_csv(CENTRALITY_FILE)
    metrics = centrality.set_index("protein").to_dict("index")

    # Use community assignments from the 0.7 confidence threshold.
    assignments = pd.read_csv(COMMUNITY_FILE)
    assignments = assignments[assignments["threshold"] == THRESHOLD]
    community_map = assignments.set_index("protein")["community"].to_dict()

    # Make the layout reproducible.
    positions = nx.spring_layout(
        graph,
        seed=42,
        weight="weight",
        k=0.8,
        iterations=100,
    )

    # Draw network edges.
    edge_x, edge_y = [], []

    for source, target in graph.edges():
        x0, y0 = positions[source]
        x1, y1 = positions[target]
        edge_x.extend([x0, x1, None])
        edge_y.extend([y0, y1, None])

    edge_trace = go.Scatter(
        x=edge_x,
        y=edge_y,
        mode="lines",
        line=dict(width=0.6, color="rgba(130,130,130,0.35)"),
        hoverinfo="skip",
        showlegend=False,
    )

    nodes = list(graph.nodes())

    node_x = [positions[node][0] for node in nodes]
    node_y = [positions[node][1] for node in nodes]

    degrees = [metrics[node]["degree"] for node in nodes]
    betweenness = [
        metrics[node]["betweenness_weighted"] for node in nodes
    ]
    communities = [community_map.get(node, 0) for node in nodes]

    hover_text = [
        (
            f"<b>{node}</b>"
            f"<br>Degree: {metrics[node]['degree']}"
            f"<br>Weighted degree: {metrics[node]['weighted_degree']:.3f}"
            f"<br>Weighted betweenness: "
            f"{metrics[node]['betweenness_weighted']:.4f}"
            f"<br>Community at threshold {THRESHOLD}: "
            f"{community_map.get(node, 'Unassigned')}"
        )
        for node in nodes
    ]

    node_trace = go.Scatter(
        x=node_x,
        y=node_y,
        mode="markers+text",
        text=nodes,
        textposition="top center",
        textfont=dict(size=8),
        hovertext=hover_text,
        hoverinfo="text",
        marker=dict(
            size=[7 + min(degree, 30) * 0.45 for degree in degrees],
            color=communities,
            colorscale="Turbo",
            showscale=True,
            colorbar=dict(title="Community"),
            line=dict(width=0.5, color="white"),
        ),
        showlegend=False,
    )

    figure = go.Figure(data=[edge_trace, node_trace])

    figure.update_layout(
        title=(
            "BioNet Explorer: Human Interferon Association Network"
            "<br><sup>STRING score ≥ 0.7 | "
            "Node size: degree | Colour: community</sup>"
        ),
        template="plotly_white",
        showlegend=False,
        hovermode="closest",
        margin=dict(l=20, r=20, t=90, b=20),
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
    )

    OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)
    figure.write_html(str(OUTPUT_FILE), include_plotlyjs=True)

    print(f"Saved: {OUTPUT_FILE}")
    print(f"Nodes visualized: {graph.number_of_nodes()}")
    print(f"Edges visualized: {graph.number_of_edges()}")
    print(f"Communities coloured: {len(set(communities))}")


if __name__ == "__main__":
    main()
