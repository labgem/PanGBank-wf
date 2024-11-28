#!/usr/bin/env python


"""Provide a command line tool to select genome from a genome try."""


import argparse
import logging
import sys
from pathlib import Path
from typing import Dict
import gzip
from treeswift import read_tree_newick
from statistics import median


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Parse genome paths and taxonomy.",
        epilog="Example: python parse_genomes_and_taxonomy.py --genomes <genomes_file> --taxonomy <taxonomy_file>",
    )
    parser.add_argument(
        "--sorted_genomes_file",
        type=Path,
        required=True,
        help="Path to a file listing genome path file sorted from best to worst ",
    )
    parser.add_argument(
        "--tree",
        type=Path,
        required=True,
        help="Newick tree of the genomes.",
    )

    parser.add_argument(
        "--number_of_genomes",
        help="Number of genomes to select from the tree.",
        default=2000,
        type=float,
    )

    parser.add_argument(
        "--disable_bar",
        type=bool,
        default=False,
        help="Disable progress bar",
    )

    parser.add_argument(
        "--cluster_indexes_file",
        help="File where cluster genome indexes are written.",
        default="cluster_indexes.txt",
        type=Path,
    )

    parser.add_argument(
        "--method",
        help="Method used to compute the distance used to cluster genomes in the tree",
        choices=["median", "max_pair"],
        default="max_pair",
        type=str,
    )

    parser.add_argument(
        "-o",
        "--output",
        help="File where filtered genomes are written.",
        default="selected_genomes_from_tree.list",
        type=Path,
    )

    parser.add_argument(
        "-l",
        "--log-level",
        help="Desired log level (default: WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="INFO",
    )

    return parser.parse_args(argv)


def write_selected_genomes_ids(selected_genomes, outfile):

    with open(outfile, "w") as fl:
        fl.write("\n".join(selected_genomes) + "\n")


def compute_leaves_median_distance(tree):
    """
    Compute for each node of the graph, the median distance between the node and its leaves
    """

    for node in tree.traverse_postorder():
        if node.is_leaf():
            node.leaf_distances = [0]
        else:
            node.leaf_distances = []
            for child in node.children:

                node.leaf_distances += [
                    d + child.edge_length for d in child.leaf_distances
                ]

        node.distance = median(node.leaf_distances)


def compute_max_leaves_pair_distance(tree):
    """
    Compute for each node of the graph, the maximal distance between the node and two of its leaves
    """
    for node in tree.traverse_postorder():
        if node.is_leaf():
            node.leaf_distances = [0]
        else:
            node.leaf_distances = []
            for child in node.children:

                node.leaf_distances += [
                    d + child.edge_length for d in child.leaf_distances
                ]

        node.distance = sum(
            sorted(node.leaf_distances, reverse=True)[
                : min(len(node.leaf_distances), 2)
            ]
        )


def cluster_tree(tree, num_cluster, method):
    """ """

    if method == "median":
        compute_leaves_median_distance(tree)
    elif method == "max_pair":
        compute_max_leaves_pair_distance(tree)
    else:
        raise ValueError(f"Unknown method to compute node distance {method}")

    sorted_nodes = sorted(
        tree.traverse_postorder(leaves=False),
        key=lambda node: (node.distance, node.num_nodes()),
        reverse=True,
    )

    node_clusters = set(tree.traverse_leaves())
    children_nodes = set()

    while len(node_clusters) > num_cluster:

        node = sorted_nodes.pop()
        if node in children_nodes:
            continue

        child_nodes = set(node.traverse_postorder())
        children_nodes |= child_nodes
        node_clusters -= child_nodes
        node_clusters.add(node)

    return node_clusters


def write_clusters(clusters, cluster_file):

    with open(cluster_file, "w") as fl:
        for indexes in clusters:
            fl.write(" ".join(map(str, sorted(indexes))) + "\n")


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)

    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")

    if not args.sorted_genomes_file.is_file():
        logging.error(f"Sorted genome file {args.sorted_genomes_file} was not found!")
        sys.exit(2)

    if not args.tree.is_file():
        logging.error(f"newick tree file {args.tree} was not found!")
        sys.exit(2)

    if not args.output.parent.exists():
        raise FileNotFoundError(
            f"Cannot write selected genomes list in '{args.output}' because its parent directory does not exists."
        )

    sorted_genomes_file = args.sorted_genomes_file

    with open(sorted_genomes_file) as fl:

        index_to_genome = {index: genome.rstrip() for index, genome in enumerate(fl)}

    tree_file = args.tree

    proper_open = gzip.open if tree_file.suffix == ".gz" else open
    with proper_open(tree_file, "rt") as fl:
        tree = read_tree_newick(fl.read().strip())

    num_cluster = args.number_of_genomes

    node_clusters = cluster_tree(tree, num_cluster, args.method)

    print(f"Found {len(node_clusters)} clusters when processing the tree.")

    selected_genomes = []
    clusters = []
    for node in node_clusters:

        cluster_indexes = [int(leaf.label) for leaf in node.traverse_leaves()]
        clusters.append(cluster_indexes)

        selected_genome_index = min(cluster_indexes)
        selected_genome = index_to_genome[selected_genome_index]
        selected_genomes.append(selected_genome)

    assert len(selected_genomes) == len(node_clusters)

    selected_genome_outfile = args.output
    print(f"Writting selected genomes in {selected_genome_outfile}.")
    write_selected_genomes_ids(selected_genomes, outfile=selected_genome_outfile)

    cluster_composition_file = args.cluster_indexes_file
    print(f"Writing clusters in {cluster_composition_file}.")
    write_clusters(clusters, cluster_composition_file)


if __name__ == "__main__":
    main()
