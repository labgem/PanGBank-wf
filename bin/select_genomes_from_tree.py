#!/usr/bin/env python

import argparse
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Set
import gzip
from treeswift import read_tree_newick, Tree
from statistics import median


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Parse genome paths and taxonomy.",
        epilog="Example: python parse_genomes_and_taxonomy.py --genomes <genomes_file> --taxonomy <taxonomy_file>",
    )
    parser.add_argument(
        "--sorted_genomes",
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
        "--cluster_composition",
        help="File where cluster genome indexes are written.",
        default="cluster_indexes.txt",
        type=Path,
    )

    parser.add_argument(
        "--selected_genomes",
        help="File where filtered genomes are written.",
        default="selected_genomes_from_tree.txt",
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
        "-l",
        "--log-level",
        help="Desired log level (default: WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="INFO",
    )

    return parser.parse_args(argv)


def write_selected_genomes(selected_genomes: List[str], outfile: Path) -> None:
    """
    Writes the selected genomes to an output file.

    :param selected_genomes: A list of selected genome paths.
    :param outfile: Path to the output file where the selected genomes will be written.
    """
    proper_open = gzip.open if outfile.suffix == ".gz" else open
    with proper_open(outfile, "wt") as fl:
        for genome_path in selected_genomes:
            fl.write(f"{genome_path}\n")


def compute_leaves_median_distance(tree: Tree) -> None:
    """
    Computes the median distance from each node to its leaves in the tree.

    :param tree: A TreeSwift `Tree` object.
    """
    for node in tree.traverse_postorder():
        if node.is_leaf():
            # Leaf nodes have a self-distance of 0
            node.leaf_distances = [0]
        else:
            # Non-leaf nodes calculate distances by summing child distances and edge lengths
            node.leaf_distances = []
            for child in node.children:
                node.leaf_distances.extend(
                    d + child.edge_length for d in child.leaf_distances
                )

        # Store the median distance of the node to its leaves
        node.distance = median(node.leaf_distances)


def compute_max_leaves_pair_distance(tree: Tree) -> None:
    """
    Computes the maximum distance between any two leaves for each node in the tree.

    :param tree: A TreeSwift `Tree` object.
    """
    for node in tree.traverse_postorder():
        if node.is_leaf():
            # Leaf nodes have a self-distance of 0
            node.leaf_distances = [0]
        else:
            # Non-leaf nodes calculate distances by summing child distances and edge lengths
            node.leaf_distances = []
            for child in node.children:
                node.leaf_distances.extend(
                    d + child.edge_length for d in child.leaf_distances
                )

        # Compute and store the maximum pairwise leaf distance
        # This is the sum of the two largest distances
        node.distance = sum(sorted(node.leaf_distances, reverse=True)[:2])


def cluster_tree(tree, num_cluster: int, method: str) -> Set:
    """
    Clusters a tree into a specified number of clusters using the chosen method
    to compute node distances.

    :param tree: The input tree to be clustered, represented in a format supporting traversal and distance computation.
    :param num_cluster: The desired number of clusters to create from the tree.
    :param method: The method used to compute node distances. Options are:
                   - "median": Uses the median distance between leaves.
                   - "max_pair": Uses the maximum pairwise distance between leaves.
    :return: A set of clustered nodes representing the final tree clusters.
    :raises ValueError: If the specified method is not recognized.
    """
    # Apply the chosen distance computation method
    if method == "median":
        compute_leaves_median_distance(tree)
    elif method == "max_pair":
        compute_max_leaves_pair_distance(tree)
    else:
        raise ValueError(f"Unknown method to compute node distance: {method}")

    # Sort nodes in descending order by distance and number of nodes
    sorted_nodes = sorted(
        tree.traverse_postorder(leaves=False),
        key=lambda node: (node.distance, node.num_nodes()),
        reverse=True,
    )

    # Initialize clusters with all leaf nodes and an empty set for processed child nodes
    node_clusters = set(tree.traverse_leaves())
    children_nodes = set()

    # Reduce the number of clusters by merging nodes
    while len(node_clusters) > num_cluster:
        node = sorted_nodes.pop()
        if node in children_nodes:
            continue

        # Merge the child nodes of the current node into the cluster
        child_nodes = set(node.traverse_postorder())
        children_nodes |= child_nodes
        node_clusters -= child_nodes
        node_clusters.add(node)

    return node_clusters


def write_clusters(clusters: List[List[int]], cluster_file: Path) -> None:
    """
    Writes the clusters to a file, with each cluster represented as a space-separated line of indexes.

    :param clusters: A list of clusters, where each cluster is a list of node indexes.
    :param cluster_file: Path to the output file where the clusters will be written.
    """
    proper_open = gzip.open if cluster_file.suffix == ".gz" else open

    with proper_open(cluster_file, "wt") as fl:
        for indexes in clusters:
            fl.write(" ".join(map(str, sorted(indexes))) + "\n")


def parse_genome_name_to_path_file(genome_name_to_path: Path) -> Dict[str, str]:
    """
    Parses a genome name-to-path file and returns a dictionary mapping paths to genome names.

    :param genome_name_to_path: Path to the file containing genome name and path mappings.
                                Each line of the file should have the format:
                                `<genome_name>\t<genome_path>`.
    :return: A dictionary where the keys are genome paths and the values are genome names.
    """
    proper_open = gzip.open if genome_name_to_path.suffix == ".gz" else open

    with proper_open(genome_name_to_path, "rt") as fl:
        path_to_genome_name = {
            line.split("\t")[1].strip(): line.split("\t")[0] for line in fl
        }

    return path_to_genome_name


def select_genome_from_tree(
    tree_file: Path, num_cluster: int, method: str, index_to_genome: Dict[int, str]
) -> Tuple[List[str], List[List[int]]]:
    """
    Selects representative genomes from a phylogenetic tree by clustering its nodes.

    :param tree_file: Path to the Newick tree file (can be gzipped or plain text).
    :param num_cluster: Number of clusters to generate from the tree.
    :param method: Clustering method to use for splitting the tree.
    :param index_to_genome: A dictionary mapping tree node indexes to genome names.
    :return: A tuple containing a list of selected genome names and a list of clusters
             (each cluster is a list of node indexes).
    """

    # Open the tree file, supporting both gzipped and plain text formats
    proper_open = gzip.open if tree_file.suffix == ".gz" else open
    with proper_open(tree_file, "rt") as fl:
        tree = read_tree_newick(fl.read().strip())

    # Cluster the tree using the specified method and number of clusters
    node_clusters = cluster_tree(tree, num_cluster, method)

    logging.info(f"Found {len(node_clusters)} clusters while processing the tree.")

    selected_genomes = []
    clusters = []

    for node in node_clusters:
        # Get the indexes of all leaves in the current cluster
        cluster_indexes = [int(leaf.label) for leaf in node.traverse_leaves()]
        clusters.append(cluster_indexes)

        # Select the genome corresponding to the smallest index in the cluster
        selected_genome_index = min(cluster_indexes)
        selected_genome = index_to_genome[selected_genome_index]
        selected_genomes.append(selected_genome)

    # Ensure the number of selected genomes matches the number of clusters
    assert len(selected_genomes) == len(node_clusters)

    return selected_genomes, clusters


def check_input_output_args(args: argparse.Namespace) -> None:
    """
    Validates the input and output file paths provided in the arguments.

    :param args: Parsed command-line arguments containing file paths for input and output.
    :raises FileNotFoundError: If any of the required input files are missing or if the output directory does not exist.
    """
    # Check if the sorted genomes file exists
    if not args.sorted_genomes.is_file():
        raise FileNotFoundError(
            f"The specified sorted genome file '{args.sorted_genomes}' does not exist."
        )

    # Check if the Newick tree file exists
    if not args.tree.is_file():
        raise FileNotFoundError(
            f"The specified Newick tree file '{args.tree}' does not exist."
        )

    # Check if the directory for the selected genomes file exists
    if not args.selected_genomes.parent.exists():
        raise FileNotFoundError(
            f"The output directory for the selected genomes file "
            f"('{args.selected_genomes.parent}') does not exist. Please create it before proceeding."
        )


def main(argv=None):
    """
    Coordinates argument parsing, program execution, and writes the selected genomes
    and clusters to specified output files.

    :param argv: Optional list of command-line arguments. Defaults to None, which uses sys.argv.
    """
    args = parse_args(argv)

    logging.basicConfig(level="INFO", format="[%(levelname)s] %(message)s")

    # Check the validity of input and output files
    check_input_output_args(args)

    # Read the sorted genomes file and prepare the genome index mapping
    sorted_genomes_file = args.sorted_genomes
    num_cluster = args.number_of_genomes
    tree_file = args.tree

    proper_open = gzip.open if sorted_genomes_file.suffix == ".gz" else open
    with proper_open(sorted_genomes_file, "rt") as fl:
        index_to_genome = {index: genome.rstrip() for index, genome in enumerate(fl)}

    # Select genomes from the tree and obtain the clusters
    selected_genomes, clusters = select_genome_from_tree(
        tree_file, num_cluster, args.method, index_to_genome
    )

    # Write selected genomes to the output file
    selected_genomes_outfile = args.selected_genomes
    logging.info(f"Writing selected genomes to {selected_genomes_outfile}.")
    write_selected_genomes(selected_genomes, outfile=selected_genomes_outfile)

    # Write the cluster composition to the output file
    logging.info(f"Writing clusters to {args.cluster_composition}.")
    write_clusters(clusters, args.cluster_composition)


if __name__ == "__main__":
    main()
