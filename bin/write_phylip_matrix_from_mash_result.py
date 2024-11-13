#!/usr/bin/env python


"""Provide a command line tool to parse genome paths and taxonomy."""


import argparse
import logging
import sys
from pathlib import Path
import gzip
import scipy.sparse
from scipy.sparse import dok_matrix
from typing import Dict

from tqdm import tqdm

from scipy.cluster.hierarchy import linkage, dendrogram, fcluster

from scipy.sparse import dok_matrix, triu, find
import numpy as np


import time
import numpy as np
import pandas as pd


def parse_mash_dist_result_into_matrix(genome_to_index:Dict[str, int], mash_result_file:Path, disable_bar:bool, multiplication_factor=10000):
    """
    """

    assert multiplication_factor <= 32767

    if not mash_result_file.is_file():
        raise FileNotFoundError(f"Mash result file '{mash_result_file}' does not exist.")

    genome_count = len(genome_to_index)

    # Using int16 to save memory (2 bytes per value)
    sparse_similarity_matrix_mash = dok_matrix((genome_count, genome_count), dtype=np.int16)

    proper_open = gzip.open if mash_result_file.suffix == ".gz" else open

    with tqdm(unit="k genome pair", disable=disable_bar, total=genome_count*genome_count / 1000) as progress:
        with proper_open(mash_result_file, "rt") as matf:
            for i, line in enumerate(matf):
                path1, path2, dist = line.split()[:3]
                index1 = genome_to_index[path1]
                index2 = genome_to_index[path2]

                # Multiply distance by 1000 and convert to integer
                similarity_int = int((1 - float(dist)) * 10000)

                if  index1 == index2:
                    pass

                elif index1 < index2:
                    sparse_similarity_matrix_mash[index1, index2] = similarity_int
                else:
                    sparse_similarity_matrix_mash[index2, index1] = similarity_int

                if i % 100000 == 0:
                    progress.update(100)

    return sparse_similarity_matrix_mash


def write_phylip_matrix(index_to_genome, sparse_similarity_matrix, phylip_matrix_file, multiplication_factor=10000):
    number_of_genome = len(index_to_genome)
    with open(phylip_matrix_file, 'w') as fl:

        fl.write(f"{number_of_genome}\n")

        for index1 in tqdm(range(number_of_genome), total=number_of_genome, unit="genome"):
            fl.write('.'.join(index_to_genome[index1].split('/')[-1].split('.')[:-2]))

            for index2 in range(index1):
                similarity = sparse_similarity_matrix[index2,index1]

                distance = 1 - (similarity / multiplication_factor)

                fl.write(f"\t{distance:.4}")
            fl.write('\n')



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
        "--mash_dist_result",
        type=Path,
        required=True,
        help="Result file of mash dist",
    )

    parser.add_argument(
        "--disable_bar",
        type=bool,
        default=False,
        help="Disable progress bar",
    )
    parser.add_argument(
        "-o",
        "--phylip_matrix",
        help="Directory where diltered genomes are stored.",
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
        raise FileNotFoundError(f"Cannot write phylip_matrix '{args.phylip_matrix}' because its parent directory does not exists.")



    multiplication_factor = 10000

    sorted_genomes_file = args.sorted_genomes_file
    with open(sorted_genomes_file) as fl:
        sorted_genomes = [line.rstrip() for line in fl]


    genome_to_index = {genome: index for index, genome in enumerate(sorted_genomes)}
    index_to_genome = {index: genome for index, genome in enumerate(sorted_genomes)}

    # sparse_matrix_file:Path = args.output / "sparse_matrix_mash_dist.npz"

    # if sparse_matrix_file.is_file():
    #     logging.info(f"Loading matrix contained in {sparse_matrix_file}")
    #     # convert matrix returned by load_npz (coo format, as saved) to dok format
    #     sparse_similarity_matrix = scipy.sparse.load_npz(sparse_matrix_file).todok()

    # else:

    logging.info(f"Parsing mash distances from '{ args.mash_dist_result}' into a sparse matrix.")

    sparse_similarity_matrix = parse_mash_dist_result_into_matrix(genome_to_index, args.mash_dist_result, disable_bar=args.disable_bar,
                                                                    multiplication_factor=multiplication_factor)

    # logging.info("Saving matrix to npz file to be loaded quicker if needed later")
    # # Convert dok_matrix to coo format, as dok format is not allowed by save_npz
    # coo_mat = sparse_similarity_matrix.tocoo()
    # scipy.sparse.save_npz(sparse_matrix_file, coo_mat)

    phylip_matrix_file = args.phylip_matrix

    logging.info(f"Writting the sparse matrix into a low triangle phylip matrix in '{phylip_matrix_file}'")
    write_phylip_matrix(index_to_genome, sparse_similarity_matrix, phylip_matrix_file, multiplication_factor=multiplication_factor)


if __name__ == "__main__":
    main()
