#!/usr/bin/env python

import argparse
import logging
from pathlib import Path
from typing import List, Dict, Tuple, Optional
import numpy as np
from collections import defaultdict
import pandas as pd


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Process clustering and distance data for analysis."
    )

    # Input arguments
    input_group = parser.add_argument_group("Input Arguments")
    input_group.add_argument(
        "--cluster_composition",
        type=Path,
        required=True,
        help="Path to a file containing cluster composition. Each line represents a cluster, with genome IDs "
        "belonging to the cluster separated by spaces.",
    )
    input_group.add_argument(
        "--phylip_matrix",
        type=Path,
        required=True,
        help="Path to the Phylip matrix file used for distance calculations.",
    )

    input_group.add_argument(
        "--species",
        type=str,
        help="Specify species name",
    )
    # Output arguments
    output_group = parser.add_argument_group("Output Arguments")
    output_group.add_argument(
        "--cluster_stat",
        type=Path,
        default="cluster_stat.tsv",
        help="Path to the output TSV file containing intra- and inter-cluster distance statistics. ",
    )
    output_group.add_argument(
        "-c",
        "--distance_count_file",
        type=Path,
        default="distance_to_count.tsv",
        help="Path to the output TSV file with counts of intra- and inter-cluster distances, "
        "used for generating density plots.",
    )

    return parser.parse_args(argv)


def count_distances(distances: List[float]):

    distances_count = defaultdict(int)

    for distance in distances:
        distances_count[np.float16(distance)] += 1

    return distances_count


def get_intra_and_inter_cluster_distances(
    phylip_matrix_file: Path, clusters: List[List[str | int]]
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculate intra- and inter-cluster distances from a Phylip matrix file.

    :param phylip_matrix_file: Path to the Phylip matrix file containing pairwise distances.
    :param clusters: List of clusters, where each cluster is a list of genome IDs.
    :return: A tuple containing two numpy arrays:
             - The first array contains intra-cluster distances.
             - The second array contains inter-cluster distances.
    """
    # Create a mapping of genome IDs to their respective cluster IDs
    genome_id_to_cluster = {
        genome_id: cluster_id
        for cluster_id, cluster in enumerate(clusters)
        for genome_id in cluster
    }
    intra_cluster_distances = []
    inter_cluster_distances = []

    # Open the Phylip matrix file and process its contents
    with open(phylip_matrix_file) as file:
        # Read the first line, which specifies the total number of genomes
        genome_count = int(file.readline().strip())
        assert len(genome_id_to_cluster) == genome_count, (
            f"Mismatch between the number of genomes in {phylip_matrix_file} "
            f"({genome_count}) and the provided clusters mapping "
            f"({len(genome_id_to_cluster)} genomes)."
        )

        # Parse the pairwise distance lines
        for i, distance_line in enumerate(file):

            distances = distance_line.strip().split()[1:]  # Skip the genome ID column

            for j, distance in enumerate(distances):
                # Classify the distance as intra- or inter-cluster
                if genome_id_to_cluster[i] == genome_id_to_cluster[j]:
                    intra_cluster_distances.append(np.float32(distance))
                else:
                    inter_cluster_distances.append(np.float32(distance))

    # Calculate expected counts for validation
    total_pairwise_count = (genome_count * (genome_count - 1)) // 2
    expected_intra_count = sum(
        (len(cluster) * (len(cluster) - 1)) // 2 for cluster in clusters
    )
    expected_inter_count = total_pairwise_count - expected_intra_count

    # Validate the calculated counts
    assert len(intra_cluster_distances) == expected_intra_count, (
        f"Expected {expected_intra_count} intra-cluster distances, but got "
        f"{len(intra_cluster_distances)}."
    )
    assert len(inter_cluster_distances) == expected_inter_count, (
        f"Expected {expected_inter_count} inter-cluster distances, but got "
        f"{len(inter_cluster_distances)}."
    )

    return np.array(intra_cluster_distances), np.array(inter_cluster_distances)


def write_distance_count_table(
    intra_distances_count: Dict[float, int],
    inter_distances_count: Dict[float, int],
    output_file: Path,
    species: str,
) -> None:
    """
    Writes a distance count table to a file, combining intra- and inter-cluster distances.

    :param intra_distances_count: A dictionary where keys are distances and values are their counts for intra-cluster distances.
    :param inter_distances_count: A dictionary where keys are distances and values are their counts for inter-cluster distances.
    :param output_file: Path to the output file where the table will be written in TSV format.
    :param species: species name to label the table
    """
    # Convert intra-cluster distances to a DataFrame
    df_intra_distances = pd.DataFrame(
        intra_distances_count.items(), columns=["distance", "count"]
    )
    df_intra_distances["type"] = "intra_cluster"

    # Convert inter-cluster distances to a DataFrame
    df_inter_distances = pd.DataFrame(
        inter_distances_count.items(), columns=["distance", "count"]
    )
    df_inter_distances["type"] = "inter_cluster"

    # Combine intra- and inter-cluster distance DataFrames
    df_combined = pd.concat([df_intra_distances, df_inter_distances])
    df_combined["species"] = species
    # Write the combined DataFrame to the output file
    df_combined.to_csv(output_file, sep="\t", index=False)


def main(argv: Optional[list[str]] = None) -> None:
    """
    Coordinate argument parsing and program execution.

    This function parses input arguments, processes a Phylip distance matrix and cluster composition file,
    computes intra- and inter-cluster distances, generates statistics, and writes results to output files.

    :param argv: List of arguments to override command-line input for testing purposes. Defaults to None.
    """
    # Parse arguments
    args = parse_args(argv)

    # Load input files
    phylip_matrix_file = args.phylip_matrix
    cluster_composition_file = args.cluster_composition

    # Read cluster composition file and construct genome-to-cluster mapping
    with open(cluster_composition_file) as fl:
        clusters = [list(map(int, line.strip().split(" "))) for line in fl if line]

    # Count single-genome clusters (singletons)
    singleton_count = len([c for c in clusters if len(c) == 1])

    # Calculate intra- and inter-cluster distances
    intra_cluster_distances, inter_cluster_distances = (
        get_intra_and_inter_cluster_distances(phylip_matrix_file, clusters)
    )

    # Count occurrences of distances
    intra_distances_count = count_distances(intra_cluster_distances)
    inter_distances_count = count_distances(inter_cluster_distances)

    # Write distance counts to file
    write_distance_count_table(
        intra_distances_count,
        inter_distances_count,
        args.distance_count_file,
        args.species,
    )

    # Compute and collect statistics
    stats = {
        "species": args.species,
        "intra_cluster_median": np.median(intra_cluster_distances),
        "intra_cluster_mean": np.mean(intra_cluster_distances),
        "intra_cluster_min": np.min(intra_cluster_distances),
        "intra_cluster_max": np.max(intra_cluster_distances),
        "intra_cluster_q1": np.percentile(intra_cluster_distances, 25),
        "intra_cluster_q3": np.percentile(intra_cluster_distances, 75),
        "inter_cluster_median": np.median(inter_cluster_distances),
        "inter_cluster_mean": np.mean(inter_cluster_distances),
        "inter_cluster_min": np.min(inter_cluster_distances),
        "inter_cluster_max": np.max(inter_cluster_distances),
        "inter_cluster_q1": np.percentile(inter_cluster_distances, 25),
        "inter_cluster_q3": np.percentile(inter_cluster_distances, 75),
        "singleton": singleton_count,
        "clusters": len(clusters),
        "genomes": sum((len(cluster) for cluster in clusters)),
    }

    # Write statistics to output file
    output_file = args.cluster_stat
    with open(output_file, "w") as out_fl:
        out_fl.write("\t".join(stats.keys()) + "\n")  # Write header
        out_fl.write("\t".join(map(str, stats.values())) + "\n")  # Write values

    logging.info(f"Cluster distance statistics saved to '{output_file}'.")


if __name__ == "__main__":
    main()
