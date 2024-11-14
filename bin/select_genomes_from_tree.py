#!/usr/bin/env python


"""Provide a command line tool to select genome from a genome try."""


import argparse
import logging
import sys
from pathlib import Path
from typing import Dict
import gzip
from treeswift import read_tree_newick


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
        help="Path to a file listing genome path file sorted from best to worst "
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
        "-o",
        "--output",
        help="Directory where filtered genomes are stored.",
        default="genome_derep_out",
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
        fl.write("\n".join(selected_genomes)+"\n")



def select_cluster_from_tree_with_max_dist(tree, num_clusters):

    # compute leaf distances and max pairwise distances
    for node in tree.traverse_postorder():
        if node.is_leaf():
            node.leaf_distance = 0
            node.max_pair_dist = 0
        else:

            node.leaf_distance = float('-inf')
            second_most_distant_leaf_dist = float('-inf')

            for child in node.children:

                current_leaf_dist = child.leaf_distance + child.edge_length

                if current_leaf_dist > node.leaf_distance:
                    second_most_distant_leaf_dist = node.leaf_distance
                    node.leaf_distance = current_leaf_dist

                elif current_leaf_dist > second_most_distant_leaf_dist:
                    second_most_distant_leaf_dist = current_leaf_dist

            node.most_distant_leaves_dist = node.leaf_distance + second_most_distant_leaf_dist

    sorted_nodes = sorted(tree.traverse_postorder(leaves=False), key=lambda node: node.most_distant_leaves_dist, reverse=True)

    node_clusters = set(tree.traverse_leaves())

    while len(node_clusters) > num_clusters:

        node = sorted_nodes.pop()

        node_clusters -= set(node.traverse_postorder())

        node_clusters.add(node)

    return node_clusters


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

    if not args.output.is_dir():
        logging.debug(f"Create output directory {args.output.absolute().as_posix()}")
        Path.mkdir(args.output, exist_ok=True)


    sorted_genomes_file = args.sorted_genomes_file
    with open(sorted_genomes_file) as fl:
        sorted_genomes = [line.rstrip() for line in fl]


    genome_to_index = {genome: index for index, genome in enumerate(sorted_genomes)}
    index_to_genome = {index: genome for index, genome in enumerate(sorted_genomes)}

    tree_file = args.tree

    proper_open = gzip.open if tree_file.suffix == ".gz" else open
    with proper_open(tree_file) as fl:
        tree = read_tree_newick(fl.read().strip())

    num_cluster = args.number_of_genomes
    node_clusters = select_cluster_from_tree_with_max_dist(tree, num_cluster)

    print(f"Found {len(node_clusters)} clusters when processing the tree.")

    selected_genomes = []
    leaf_is_index = True

    for node in node_clusters:

        if leaf_is_index:
            selected_genome_index = min([int(leaf.label) for leaf in node.traverse_leaves()])

            selected_genome = index_to_genome[selected_genome_index]
            selected_genomes.append(selected_genome)

        else:
            genomes = [leaf.label for leaf in node.traverse_leaves()]

            try:
                selected_genome = min(genomes, key=lambda x: genome_to_index[x])
            except KeyError as err:
                raise KeyError(f"Genome {err} not found in the sorted genome list file '{args.sorted_genomes_file}'. ")

            selected_genomes.append(selected_genome)
    assert len(selected_genomes) == len(node_clusters)

    selected_genome_outfile = args.output / "selected_genomes_from_tree.list"
    write_selected_genomes_ids(selected_genomes, outfile=selected_genome_outfile)


if __name__ == "__main__":
    main()
