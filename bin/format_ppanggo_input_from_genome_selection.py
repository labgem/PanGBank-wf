#!/usr/bin/env python

import argparse
import logging
from pathlib import Path
from typing import Dict


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Parse genome paths and taxonomy.",
        epilog="Example: python parse_genomes_and_taxonomy.py --genomes <genomes_file> --taxonomy <taxonomy_file>",
    )

    parser.add_argument(
        "--selected_genomes",
        help="Path to a file containing a list of selected genome paths, one per line. ",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--genome_name_to_path",
        type=Path,
        required=True,
        help="Path to a TSV file mapping genome names to genome file paths. "
        "The file must have two columns: the first column for genome names, "
        "and the second column for the corresponding genome file paths.",
    )

    parser.add_argument(
        "--fasta_to_original_path",
        type=Path,
        required=False,
        help="Path to a TSV file mapping temporary FASTA file paths to their original genome file paths. "
        "This file must have two columns: the first column for temporary FASTA file paths "
        "found in the --genome_name_to_path and --selected_genomes files, "
        "and the second column for the original genome file paths to be used in the output.",
    )

    parser.add_argument(
        "--formatted_genomes",
        type=Path,
        default="ppanggo_input_selected_genomes.tsv",
        help="Path to a TSV file mapping temporary FASTA file paths to their original genome file paths. "
        "This file must have two columns: the first column for temporary FASTA file paths "
        "found in the --genome_name_to_path and --selected_genomes files, "
        "and the second column for the original genome file paths to be used in the output.",
    )

    return parser.parse_args(argv)


def parse_genome_name_to_path_file(genome_name_to_path: Path) -> Dict[str, str]:
    """
    Parses a genome name-to-path file and returns a dictionary mapping paths to genome names.

    :param genome_name_to_path: Path to the file containing genome name and path mappings.
                                Each line of the file should have the format:
                                `<genome_name>\t<genome_path>`.
    :return: A dictionary where the keys are genome paths and the values are genome names.
    """
    with open(genome_name_to_path, "r") as fl:
        path_to_genome_name = {
            line.split("\t")[1].strip(): line.split("\t")[0] for line in fl
        }

    return path_to_genome_name


def main(argv=None):
    """
    Coordinates argument parsing, program execution, and writes the selected genomes
    and clusters to specified output files.

    :param argv: Optional list of command-line arguments. Defaults to None, which uses sys.argv.
    """
    args = parse_args(argv)

    logging.basicConfig(level="INFO", format="[%(levelname)s] %(message)s")

    # Read the sorted genomes file and prepare the genome index mapping
    selected_genomes_file = args.selected_genomes

    fasta_to_original_path = args.fasta_to_original_path

    with open(selected_genomes_file) as fl:
        selected_genomes = [genome_path.rstrip() for genome_path in fl]

    tmp_fasta_to_original_path = {}
    if args.fasta_to_original_path:
        with open(fasta_to_original_path) as fl:
            tmp_fasta_to_original_path = {
                line.split("\t")[0]: line.split("\t")[1].rstrip() for line in fl
            }

    # Parse the genome name to path mapping
    path_to_genome_name = parse_genome_name_to_path_file(args.genome_name_to_path)

    # Write selected genomes to the output file
    selected_genomes_outfile = args.selected_genomes
    logging.info(f"Writing selected genomes to {selected_genomes_outfile}.")

    with open(args.formatted_genomes, "w") as fl:
        for genome_path in selected_genomes:

            genome_name = path_to_genome_name[genome_path]

            if tmp_fasta_to_original_path:
                genome_path = tmp_fasta_to_original_path[genome_path]

            fl.write(f"{genome_name}\t{genome_path}\n")


if __name__ == "__main__":
    main()
