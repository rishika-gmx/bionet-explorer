# BioNet Explorer 

**Network analysis and visualization of human interferon-associated protein interactions**

BioNet Explorer is a computational biology project that explores the organization of human interferon-associated protein networks using public biological resources, graph theory, and interactive visualization.

The project combines pathway-derived gene curation, STRING functional association data, network centrality analysis, community detection, and confidence-threshold sensitivity analysis to investigate how interferon signaling proteins connect with antiviral defense and other cellular processes.

## Research Question

**How is the human interferon-associated network organized, and which proteins occupy highly connected or bridging positions between functional groups?**

## Objectives

- Curate human protein candidates from Reactome's Human Interferon Signaling pathway participants.
- Map candidate genes to STRING protein identifiers and retrieve functional associations.
- Construct and validate a protein association network using Python and NetworkX.
- Identify highly connected proteins and potential network bridges using centrality measures.
- Detect functional communities using Louvain community detection.
- Assess the sensitivity of network structure and community membership to different STRING confidence thresholds.
- Generate interactive visualizations of the complete network and a stable interferon/antiviral-associated module.

## Tools and Data Sources

| Tool or resource | Purpose |
|---|---|
| Python | Data processing and analysis |
| pandas | Data cleaning, tabulation, and result management |
| NetworkX | Graph construction, centrality, and community detection |
| Plotly | Interactive network visualization |
| STRING | Protein functional association data and confidence scores |
| Reactome | Human interferon-signaling pathway participants |
| UniProt | Protein reference and species validation |
| Git and GitHub | Version control and project sharing |

## Workflow

1. **Pathway-derived gene curation:** Retrieve Reactome pathway participants and identify human protein references using UniProt species information.
2. **Identifier mapping:** Map curated gene symbols to STRING protein identifiers and retain the unmapped candidates in the curation records.
3. **Network construction:** Retrieve STRING associations at a minimum confidence score of 0.700 and construct an undirected graph.
4. **Validation:** Check interaction records for missing endpoints, invalid scores, duplicate edges, self-loops, and nonhuman records.
5. **Centrality analysis:** Calculate degree, weighted degree, and betweenness centrality using STRING confidence scores to define edge distances.
6. **Community detection:** Apply the Louvain algorithm at confidence thresholds of 0.700 and 0.900.
7. **Sensitivity analysis:** Compare community partitions using the Adjusted Rand Index (ARI) and compare selected module memberships using Jaccard similarity.
8. **Visualization:** Generate interactive HTML visualizations of the full network and the shared interferon/antiviral-associated module.

## Key Results

### Network overview

| Metric | Result |
|---|---:|
| Human candidate gene symbols | 338 |
| Uniquely STRING-mapped candidates | 335 |
| Unmapped candidates | 3 |
| Interaction records | 4,486 |
| Proteins in the resulting network | 324 |
| Connected components | 2 |
| Network density | 0.086 |
| STRING confidence threshold | ≥ 0.700 |

The interaction table passed the implemented structural validation checks, with no missing endpoint rows, invalid scores, nonhuman rows, duplicate edges, or self-loops reported.

### Centrality analysis

STAT1 had the highest weighted betweenness centrality in the resulting network.

| Protein | Degree | Weighted betweenness |
|---|---:|---:|
| STAT1 | 108 | 0.1842 |
| TP53 | 63 | 0.1054 |
| MAPT | 86 | 0.1002 |
| EIF2AK2 | 76 | 0.0987 |
| RPS27A | 92 | 0.0747 |

These results identify proteins with prominent topological positions. High centrality alone does not establish biological importance, causal regulation, or pathway activity.

### Community detection and stability

| Metric | Threshold 0.700 | Threshold 0.900 |
|---|---:|---:|
| Proteins | 324 | 306 |
| Communities detected | 7 | 10 |
| Modularity | 0.6003 | 0.6135 |

The Adjusted Rand Index between the two full network partitions was **0.6607**, calculated across 306 proteins shared between the partitions.

An interferon/antiviral-associated community was identified at both thresholds. The two communities shared 74 proteins, with a Jaccard similarity of **0.6667**. The induced subnetwork of these shared proteins contained 913 edges and formed one connected component.

Shared members included STAT1, STAT2, JAK1, JAK2, TYK2, IRF9, ISG15, DDX58, MAVS, MX1, and members of the OAS and IFIT families.

Together, these findings show that the inferred network contains recognizable interferon-associated functional structure, while also demonstrating that community assignments depend partly on the selected confidence threshold.

## Interactive Visualizations

The project generates two interactive HTML files:

- **Full network:** `results/visualizations/reference_network/interferon_network.html`
- **Stable interferon/antiviral-associated module:** `results/visualizations/reference_network/stable_interferon_module.html`

Open either file in a web browser to explore protein labels, connections, and network metrics.

## Repository Structure

```text
bionet-explorer/
├── data/
│   ├── interferon_genes.txt
│   └── reference/
├── src/
│   ├── string_api.py
│   ├── validation.py
│   ├── network_analysis.py
│   ├── community_analysis.py
│   ├── community_stability.py
│   ├── sensitivity_analysis.py
│   ├── curate_reactome_genes.py
│   ├── fetch_curated_network.py
│   ├── visualize_curated_network.py
│   └── visualize_stable_module.py
├── results/
├── visualizations/
│   ├── interferon_network.html
│   └── reference_network/
│       ├── interferon_network.html
│       └── stable_interferon_module.html
│   └── reference_network/
├── tests/
├── app.py
├── requirements.txt
├── .gitignore
└── README.md
```

## Running the Project

### Requirements

- Python 3
- Git
- Internet access for retrieving data from external APIs

### Setup

Clone the repository and enter the project directory:

```bash
git clone <YOUR-GITHUB-REPOSITORY-URL>
cd bionet-explorer
```

Create and activate a virtual environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
```

Install the dependencies:

```bash
pip install -r requirements.txt
```

### Generate the curated network visualizations

The commands below assume that the curated interaction data and analysis outputs are already present.

Run the full-network visualization:

```bash
python src/visualize_curated_network.py
```

Run the focused module visualization:

```bash
python src/visualize_stable_module.py
```

The generated HTML files are saved in `results/reference_network/`.

For a fresh run, execute the data retrieval, curation, validation, and analysis stages in the appropriate order. API access and response contents may change over time.

## Limitations

- STRING associations represent functional associations and are not necessarily direct physical protein–protein interactions.
- The Reactome-derived candidate set is based on extracted pathway participants and is not guaranteed to represent an exhaustive list of core interferon-signaling genes.
- The resulting network contains fewer distinct proteins than the number of uniquely mapped input candidates; this difference should be accounted for when interpreting network coverage.
- Centrality and community structure depend on network construction choices, confidence thresholds, and algorithm settings.
- The stable module is defined computationally from shared community membership at two thresholds. It has not been independently validated experimentally.
- Network topology alone cannot establish gene expression, protein activity, causal relationships, or clinical significance.

## Future Directions

- Improve visualization readability through focused network views and selective node labeling.
- Investigate the biological context of high-centrality proteins and detected communities.
- Improve reproducibility with documented analysis commands and automated tests.
- Explore an interactive interface for inspecting network structure and protein-level metrics.

## Acknowledgements

This project uses publicly available biological resources from Reactome, UniProt, and STRING.

---

*BioNet Explorer is an educational computational biology portfolio project exploring biological network analysis and visualization.*
