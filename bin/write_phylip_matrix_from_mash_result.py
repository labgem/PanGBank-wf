#!/usr/bin/env python


"""Provide a command line tool to parse genome paths and taxonomy."""


import argparse
import logging
import sys
from pathlib import Path
import gzip

# from typing import Dict

from tqdm import tqdm

import numpy as np
from collections import defaultdict
import pandas as pd


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
        "--mash_dist_result",
        type=Path,
        required=True,
        help="Result file of mash dist",
    )

    parser.add_argument(
        "--disable_bar",
        action="store_true",
        default=False,
        help="Disable progress bar",
    )
    parser.add_argument(
        "-o",
        "--phylip_matrix",
        help="phylip matrix file.",
        default="distance.phylip",
        type=Path,
    )
    parser.add_argument(
        "-c",
        "--distance_count_file",
        help="TSV file with count of ditance.",
        default="distance_to_count.tsv",
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


def initialise_lower_triangular_list_of_np_array(size):

    list_of_array = [np.ones(row, dtype=np.float32) for row in range(size)]

    return list_of_array


def write_distance_count_table(distances_count, output_file):

    df_count = pd.DataFrame(distances_count.items(), columns=["dist", "count"])
    df_count.to_csv(output_file, sep="\t", index=False)


def parse_mash_dist_result_np_array(
    genome_to_index, mash_result_file, disable_bar=False
):

    distances_count = defaultdict(int)
    genome_count = len(genome_to_index)
    list_of_arrays = initialise_lower_triangular_list_of_np_array(genome_count)

    proper_open = gzip.open if mash_result_file.suffix == ".gz" else open

    with tqdm(
        unit="k genome pair",
        disable=disable_bar,
        total=genome_count * genome_count / 1000,
    ) as progress:
        with proper_open(mash_result_file, "rt") as matf:
            for i, line in enumerate(matf):
                path1, path2, dist = line.split()[:3]
                index1 = genome_to_index[path1]
                index2 = genome_to_index[path2]

                if index1 < index2:
                    list_of_arrays[index2][index1] = np.float32(dist)

                    distances_count[np.float16(dist)] += 1

                if i % 100000 == 0:
                    progress.update(100)

    return list_of_arrays, distances_count


def write_phylip_matrix_from_list_of_arrays(
    number_of_genomes, list_of_arrays, phylip_matrix_file, disable_bar=False
):

    with open(phylip_matrix_file, "w") as fl:
        # Write the number of genomes
        fl.write(f"{number_of_genomes}\n0\n")

        for index1 in tqdm(
            range(1, number_of_genomes),
            total=number_of_genomes,
            unit="genome",
            disable=disable_bar,
        ):

            dist_line = "\t".join((f"{d:.8f}" for d in list_of_arrays[index1]))

            fl.write(f"{index1}\t{dist_line}\n")


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)

    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")

    if not args.sorted_genomes_file.is_file():
        logging.error(f"Sorted genome file {args.sorted_genomes_file} was not found!")
        sys.exit(2)

    if not args.mash_dist_result.is_file():
        logging.error(f"mash_dist_result file {args.mash_dist_result} was not found!")
        sys.exit(2)

    if not args.phylip_matrix.parent.exists():
        raise FileNotFoundError(
            f"Cannot write phylip_matrix '{args.phylip_matrix}' because its parent directory does not exists."
        )

    if not args.distance_count_file.parent.exists():
        raise FileNotFoundError(
            f"Cannot write distance count in '{args.distance_count_file}' because its parent directory does not exists."
        )

    sorted_genomes_file = args.sorted_genomes_file
    phylip_matrix_file = args.phylip_matrix

    logging.info(f"Parsing sorted genomes from {args.sorted_genomes_file}.")
    with open(sorted_genomes_file) as fl:
        genome_to_index = {genome.rstrip(): index for index, genome in enumerate(fl)}

    logging.info(
        f"Parsing mash distances from '{ args.mash_dist_result}' into a numpy flat array."
    )
    list_of_arrays, distances_count = parse_mash_dist_result_np_array(
        genome_to_index, args.mash_dist_result, disable_bar=args.disable_bar
    )

    write_distance_count_table(distances_count, args.distance_count_file)

    logging.info(
        f"Converting the numpy flat array into a low triangle phylip matrix in '{phylip_matrix_file}'"
    )
    write_phylip_matrix_from_list_of_arrays(
        len(genome_to_index), list_of_arrays, phylip_matrix_file
    )


if __name__ == "__main__":
    main()
