
import networkx as nx
import pandas as pd
import plotly.graph_objects as go

from network_analysis import build_network

RESULTS_DIR = "results"


def main():
    graph = build_network()

    centrality = pd.read_csv(
        f"{RESULTS_DIR}/centrality_results.csv"
    )
    metrics = centrality.set_index("protein").to_dict("index")

    # Fixed seed makes the layout reproducible.
    positions = nx.spring_layout(
        graph,
        seed=42,
        weight="weight",
        k=0.8,
        iterations=100,
    )

    # Draw all edges in one trace.
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
        line=dict(width=0.7, color="rgba(130,130,130,0.45)"),
        hoverinfo="skip",
        showlegend=False,
    )

    nodes = list(graph.nodes())

    node_x = [positions[node][0] for node in nodes]
    node_y = [positions[node][1] for node in nodes]
    degrees = [metrics[node]["degree"] for node in nodes]
    betweenness = [
    metrics[node]["betweenness_unweighted"]
    for node in nodes
]

    hover_text = [
        (
            f"<b>{node}</b>"
            f"<br>Degree: {metrics[node]['degree']}"
            f"<br>Weighted degree: "
            f"{metrics[node]['weighted_degree']:.3f}"
            f"<br>Betweenness (unweighted): "
            f"{metrics[node]['betweenness_unweighted']:.4f}"
        )
        for node in nodes
    ]

    node_trace = go.Scatter(
        x=node_x,
        y=node_y,
        mode="markers+text",
        text=nodes,
        textposition="top center",
        textfont=dict(size=10),
        hovertext=hover_text,
        hoverinfo="text",
        marker=dict(
            size=[12 + degree * 1.2 for degree in degrees],
            color=betweenness,
            colorscale="Viridis",
            showscale=True,
            colorbar=dict(title="Betweenness"),
            line=dict(width=1, color="white"),
        ),
        showlegend=False,
    )

    figure = go.Figure(data=[edge_trace, node_trace])

    figure.update_layout(
        title="BioNet Explorer: Interferon Association Network",
        template="plotly_white",
        showlegend=False,
        hovermode="closest",
        margin=dict(l=20, r=20, t=70, b=20),
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
    )

    figure.write_html(
        f"{RESULTS_DIR}/visualizations/interferon_network.html",
        include_plotlyjs=True,
    )

    print("Saved results/visualizations/interferon_network.html")


if __name__ == "__main__":
    main()
