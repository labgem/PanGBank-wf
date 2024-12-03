#!/usr/bin/env python

import argparse
import logging
import sys
from pathlib import Path
import gzip
from typing import List, Dict
from tqdm import tqdm
import numpy as np


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
        type=lambda x: sys.stdin if x == "-" else Path(x),
        required=True,
        help=(
            "Result of the `mash dist` command. Provide a file path or use '-' to read from "
            "standard input (e.g., when piping the output directly)."
        ),
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


def initialise_lower_triangular_list_of_np_array(size: int) -> list[np.ndarray]:
    """
    Initializes a lower triangular list of numpy arrays, where each array contains `1`s.

    :param size: The size of the list and the maximum index of the arrays.
    :return: A list of numpy arrays, where the i-th array has `i` elements initialized to 1.
    """
    return [np.ones(row, dtype=np.float32) for row in range(size)]


def parse_mash_dist_result_np_array(
    genome_to_index: Dict[str, int], mash_dist_result: Path, disable_bar: bool = False
) -> list[np.ndarray]:
    """
    Parses a Mash distance result file and fills a lower triangular numpy array with distances.

    :param genome_to_index: A dictionary mapping genome paths to their corresponding indices.
    :param mash_dist_result: The path to the Mash distance result file (can be gzipped).
    :param disable_bar: A flag to disable the progress bar. Defaults to False.
    :return: A list of numpy arrays representing the lower triangular distance matrix.
    """
    genome_count = len(genome_to_index)
    list_of_arrays = initialise_lower_triangular_list_of_np_array(genome_count)

    # Determine the input stream based on file type
    if mash_dist_result == sys.stdin:
        mash_dist_source = sys.stdin
    else:
        proper_open = gzip.open if mash_dist_result.suffix == ".gz" else open
        mash_dist_source = proper_open(mash_dist_result, "rt")

    # Progress bar for processing pairs of genomes
    with tqdm(
        unit="k genome pair",
        disable=disable_bar,
        total=(genome_count * genome_count) / 1000,
    ) as progress:
        with mash_dist_source as input_stream:
            for i, line in enumerate(input_stream):
                path1, path2, dist = line.split()[:3]
                index1 = genome_to_index[path1]
                index2 = genome_to_index[path2]

                # Only store the lower triangular part of the distance matrix
                if index1 < index2:
                    list_of_arrays[index2][index1] = np.float32(dist)

                # Update progress bar every 100,000 lines
                if i % 100000 == 0:
                    progress.update(100)

    return list_of_arrays


def write_phylip_matrix_from_list_of_arrays(
    number_of_genomes: int,
    list_of_arrays: List[List[float]],
    phylip_matrix_file: str,
    disable_bar: bool = False,
) -> None:
    """
    Writes a Phylip distance matrix to a file from a list of arrays.

    :param number_of_genomes: The number of genomes (rows in the matrix).
    :param list_of_arrays: A list of arrays where each array represents the distance of a genome
                            to all others. Each array corresponds to a genome.
    :param phylip_matrix_file: The file path where the Phylip matrix will be written.
    :param disable_bar: A flag to disable the progress bar.
    """
    with open(phylip_matrix_file, "w") as fl:
        # Write the number of genomes and an initial '0' (indicating no species for the header)
        fl.write(f"{number_of_genomes}\n0\n")

        # Write the distance matrix, excluding the first genome as it's used for indexing
        for index1 in tqdm(
            range(1, number_of_genomes),
            total=number_of_genomes,
            unit="genome",
            disable=disable_bar,
        ):
            dist_line = "\t".join(f"{d:.8f}" for d in list_of_arrays[index1])
            fl.write(f"{index1}\t{dist_line}\n")


def check_input_output_args(args: argparse.Namespace) -> None:
    """
    Validates the input and output file paths provided in the arguments.

    :param args: Parsed command-line arguments containing file paths for input and output.
    :raises FileNotFoundError: If any of the required input files are missing or if the output directory does not exist.
    """
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")

    if not args.sorted_genomes_file.is_file():
        raise FileNotFoundError(
            f"Sorted genome file {args.sorted_genomes_file} was not found!"
        )

    if args.mash_dist_result != sys.stdin and not args.mash_dist_result.is_file():
        raise FileNotFoundError(
            f"mash_dist_result file {args.mash_dist_result} was not found!"
        )

    if not args.phylip_matrix.parent.exists():
        raise FileNotFoundError(
            f"Cannot write phylip_matrix '{args.phylip_matrix}' because its parent directory does not exists."
        )

    if not args.distance_count_file.parent.exists():
        raise FileNotFoundError(
            f"Cannot write distance count in '{args.distance_count_file}' because its parent directory does not exists."
        )


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    check_input_output_args(args)

    sorted_genomes_file = args.sorted_genomes_file
    phylip_matrix_file = args.phylip_matrix

    logging.info(f"Parsing sorted genomes from {args.sorted_genomes_file}.")
    with open(sorted_genomes_file) as fl:
        genome_to_index = {genome.rstrip(): index for index, genome in enumerate(fl)}

    logging.info(
        f"Parsing mash distances from '{ args.mash_dist_result}' into a numpy flat array."
    )
    list_of_arrays = parse_mash_dist_result_np_array(
        genome_to_index, args.mash_dist_result, disable_bar=args.disable_bar
    )

    logging.info(
        f"Converting the numpy flat array into a low triangle phylip matrix in '{phylip_matrix_file}'"
    )

    write_phylip_matrix_from_list_of_arrays(
        len(genome_to_index), list_of_arrays, phylip_matrix_file
    )


if __name__ == "__main__":
    main()
