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

def remove_close_genomes(genomes_under_selection, genome_to_index, mash_dist_matrix, genomes_removed, identity_cutoff):
    """
    """
    # Get last element (which is the 'best' genome), and remove it from the list
    ref_name = genomes_under_selection.pop()

    ref_num = genome_to_index[ref_name]

    others = genomes_under_selection[::-1]

    # For each genome, compare its distance to reference genome 'ref_name'
    for gname in others:
        # Column of genome in mat_sp
        other_num = genome_to_index[gname]
        # Get distance between the 2 genomes
        if  ref_num < other_num:
            ident = mash_dist_matrix[ref_num, other_num]
        else:
            print("Should never happen as mat_sp is a triangle matrix!")
            ident = mash_dist_matrix[other_num, ref_num]

        # If distance not in the limits, remove genome from to_try and add to genomes_removed list
        if ident > identity_cutoff:
            genomes_under_selection.remove(gname)
            genomes_removed[gname] = [ref_name, identity_cutoff]

    return 0


def parse_mash_dist_result_into_matrix(genome_to_index:Dict[str, int], mash_result_file:Path, disable_bar:bool):
    """
    """

    if not mash_result_file.is_file():
        print(f"Matrix file {mash_result_file} does not exist. We cannot read it "
                     "and do the next steps. Program ending.")
        sys.exit(1)

    genome_count = len(genome_to_index)

    # Create square matrix with nbgen cols/lines. dok format is a 'Dictionary Of Keys'
    # -> writes (0, 1) value
    sparse_matrix_mash = dok_matrix((genome_count, genome_count), dtype=float)
    # Write matrix values
    proper_open = gzip.open if mash_result_file.suffix == ".gz" else open

    with tqdm(unit="k genome pair", disable=disable_bar) as progress:
        with proper_open(mash_result_file, "rt") as matf:
            for i, line in enumerate(matf):
                path1, path2, dist = line.split()[:3]
                num1 = genome_to_index[path1]
                num2 = genome_to_index[path2]

                # only in lower triangle (no duplicate)
                if  num1 == num2:
                    pass
                elif num1 < num2:
                    sparse_matrix_mash[num1, num2] = 1 - float(dist)
                else:
                    sparse_matrix_mash[num2, num1] = 1 - float(dist)

                if i % 10000 == 0:
                    progress.update(10)

    return sparse_matrix_mash

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
        "--min_dist",
        help="Discard genome(s) closer than a Mash distance. ",
        default=0.005,
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


    identity_cutoff = 1 - args.min_dist

    sorted_genomes_file = args.sorted_genomes_file
    with open(sorted_genomes_file) as fl:
        sorted_genomes = [line.rstrip() for line in fl]


    genome_to_index = {genome: index for index, genome in enumerate(sorted_genomes)}

    sparse_matrix_file:Path = args.output / "sparse_matrix_mash_dist.npz"

    if sparse_matrix_file.is_file():
        logging.info(f"Loading matrix contained in {sparse_matrix_file}")
        # convert matrix returned by load_npz (coo format, as saved) to dok format
        mash_dist_matrix = scipy.sparse.load_npz(sparse_matrix_file).todok()

    else:
        mash_dist_matrix = parse_mash_dist_result_into_matrix(genome_to_index, args.mash_dist_result, disable_bar=args.disable_bar)

        logging.info("Saving matrix to npz file to be loaded quicker if needed later")
        # Convert dok_matrix to coo format, as dok format is not allowed by save_npz
        coo_mat = mash_dist_matrix.tocoo()
        scipy.sparse.save_npz(sparse_matrix_file, coo_mat)


    genomes_under_selection = sorted_genomes[::-1]
    # Put list of genomes removed by mash comparison, and why
    # (out of limits distance with which genome)
    genomes_removed = {}  # {genome: [compared_with, dist]}
    genome_count = len(genomes_under_selection)

    last_genome_process_count = 0

    with tqdm(total=genome_count, unit="genome", disable=args.disable_bar) as progress:

        while len(genomes_under_selection) > 1:

            remove_close_genomes(genomes_under_selection, genome_to_index, mash_dist_matrix, genomes_removed, identity_cutoff)

            current_genome_process_count = genome_count - len(genomes_under_selection)
            genome_process_in_loop = current_genome_process_count - last_genome_process_count
            last_genome_process_count = current_genome_process_count
            progress.update(genome_process_in_loop)



    logging.info("Final number of genomes in dataset: {}".format(genome_count - len(genomes_removed)))

    selected_genomes = [genome for genome in sorted_genomes_file if genome not in genomes_removed]
    selected_genome_outfile = args.output / "selected_genomes.list"

    with open(selected_genome_outfile, "w") as fl:
        fl.write("\n".join(selected_genomes))


    return genomes_removed



if __name__ == "__main__":
    main()
