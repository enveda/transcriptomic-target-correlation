# -*- coding: utf-8 -*-

"""Visualize overlap."""

import logging
from typing import List, Dict

import igraph
import pandas as pd

logger = logging.getLogger(__name__)

COLORS_CATEGORY = {
    'present_in_1_but_not_in_2': "#BEBEBE",  # grey
    'matching': "#65c281",  # green
    'non_matching': "#faa39b",  # red
    'target': "#000000",  # black
    'present_in_2_but_not_in_1': "#BEBEBE",  # grey
}

BASE_VISUAL_STYLE = {
    "bbox": (3000, 3000),
    "margin": 300,
    "vertex_size": 18,
    "vertex_label_size": 20,
    "edge_curved": False,
}


def render_bel_with_igraph(
    path: str,
    node_categories: Dict[str, str],
    gene_names: List[str],
    df: pd.DataFrame,
):
    """Plot overlap."""
    nodes_in_network = set(df.source.unique()).union(set(df.target.unique()))

    # Dict mapping each node to its color (see COLORS_CATEGORY variable above)
    node_colors = [
        COLORS_CATEGORY[node_categories[gene]]
        for gene in gene_names
        if gene in nodes_in_network
    ]

    # List of all protein HGNC symbols

    graph: igraph.Graph = igraph.Graph()
    graph.add_vertices([
        gene
        for gene in gene_names
        if gene in nodes_in_network
    ])
    graph.vs["color"] = node_colors
    # graph.vs["label"] = names

    # Add edges to the graph
    for source, target, _, _ in df.values:
        graph.add_edge(
            source,
            target,
            weight=1
        )
        graph.add_edge(
            target,
            source,
            weight=1
        )

    visual_style = BASE_VISUAL_STYLE.copy()

    # Get largest component
    # largest = graph.clusters().giant()
    # print(graph.vcount())
    # print(largest.vcount())
    visual_style["layout"] = graph.layout('fruchterman_reingold', niter=2000)

    # Plot the graph
    igraph.plot(graph, path, **visual_style)

# if __name__ == '__main__':
#     logging.basicConfig(level=logging.INFO)
#     render_bel_with_igraph(path='MY EXPORT PATH')
