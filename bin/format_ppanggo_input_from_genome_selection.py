#!/usr/bin/env python

import argparse
import logging
from pathlib import Path
from typing import Dict, Set
import gzip

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
        "--reference_genomes",
        type=Path,
        required=False,
        help="Path to a file containing a list of reference genome names, one per line. "
        "These genomes will be included in the selection regardless if they are in the provided selection.",
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
        default="ppanggo_input_selected_genomes.tsv.gz",
        help="Path to a TSV file mapping temporary FASTA file paths to their original genome file paths. "
        "This file must have two columns: the first column for temporary FASTA file paths "
        "found in the --genome_name_to_path and --selected_genomes files, "
        "and the second column for the original genome file paths to be used in the output.",
    )

    parser.add_argument(
        "--summary_selection",
        help="File where genome selection summary stats are written.",
        default="summary_selection.yaml",
        type=Path,
    )

    parser.add_argument(
        "--species",
        type=str,
        help="Specify species name",
        default=None
    )

    return parser.parse_args(argv)

def parse_genome_name_to_path_file(genome_name_to_path_file: Path) -> Dict[str, str]:
    """
    Parses a genome name-to-path file and returns a dictionary mapping paths to genome names.

    :param genome_name_to_path_file: Path to the file containing genome name and path mappings.
                                Each line of the file should have the format:
                                `<genome_name>\t<genome_path>`.
    :return: A dictionary where the keys are genome paths and the values are genome names.
    """

    proper_open = gzip.open if genome_name_to_path_file.suffix == ".gz" else open
    with proper_open(genome_name_to_path_file, "rt") as fl:
        path_to_genome_name = {
            line.split("\t")[1].strip(): line.split("\t")[0] for line in fl
        }

    return path_to_genome_name

def get_path_of_reference_genomes(
    reference_genome_names: Set[str],
    path_to_genome_name: Dict[str, str]
) -> Set[str]:
    """
    Retrieve the file paths of reference genomes that match the specified genome names.

    :param reference_genome_names: A set of genome names to be used as reference.
    :param path_to_genome_name: A mapping of file paths to genome names.
    :return: A set of file paths corresponding to the reference genomes.
    """
    # Identify reference genome names present in the provided genome name mapping
    current_species_ref_genomes = set(path_to_genome_name.values()) & reference_genome_names

    # Collect file paths corresponding to the identified reference genome names
    reference_genome_paths = {path for path, name in path_to_genome_name.items() if name in current_species_ref_genomes}

    return reference_genome_paths


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

    proper_open = gzip.open if selected_genomes_file.suffix == ".gz" else open

    with proper_open(selected_genomes_file, "rt") as fl:
        selected_genomes = {genome_path.rstrip() for genome_path in fl}

    tmp_fasta_to_original_path = {}
    if args.fasta_to_original_path:
        proper_open = gzip.open if fasta_to_original_path.suffix == ".gz" else open
        with proper_open(fasta_to_original_path, "rt") as fl:
            tmp_fasta_to_original_path = {
                line.split("\t")[0]: line.split("\t")[1].rstrip() for line in fl
            }

    # Parse the genome name to path mapping
    path_to_genome_name = parse_genome_name_to_path_file(args.genome_name_to_path)

    reference_genome_count = 0
    if args.reference_genomes:
        proper_open = gzip.open if args.reference_genomes.suffix == ".gz" else open

        with proper_open(args.reference_genomes, "rt") as fl:
            # Load reference genome names from the file
            reference_genomes = {line.rstrip() for line in fl if line}
            logging.info(f"Loaded {len(reference_genomes)} reference genome names from {args.reference_genomes}.")

        # Get paths for the specified reference genomes
        reference_genome_paths = get_path_of_reference_genomes(reference_genomes, path_to_genome_name)
        logging.info(f"Retrieved {len(reference_genome_paths)} reference genome paths matching the provided input genomes.")

        # Add the retrieved genome paths to the selected genomes set
        selection_count = len(selected_genomes)
        reference_genome_count = len(reference_genome_paths)
        selected_genomes |= reference_genome_paths
        new_selection_count = len(selected_genomes)
        added_genomes_count = new_selection_count - selection_count

        logging.info(
            f"Added {added_genomes_count} new reference genomes to the selected list, "
            f"bringing the total to {new_selection_count}. These genomes were not already in the list."
        )

    # Write selected genomes to the output file
    logging.info(f"Writing selected genomes to {args.formatted_genomes}")

    proper_open = gzip.open if args.formatted_genomes.suffix == ".gz" else open

    with proper_open(args.formatted_genomes, "wt") as fl:
        for genome_path in selected_genomes:
            genome_name = path_to_genome_name[genome_path]

            if tmp_fasta_to_original_path:
                genome_path = tmp_fasta_to_original_path[genome_path]

            fl.write(f"{genome_name}\t{genome_path}\n")

    # Write the cluster composition to the output file
    logging.info(f"Writing summary selection to {args.summary_selection}.")

    final_selection = len(selected_genomes)
    selection_no_ref =  final_selection - reference_genome_count
    total_genomes = len(path_to_genome_name)
    discarded_genome_count = total_genomes - final_selection

    with open(args.summary_selection, "w") as fl:

        fl.write('species\tselected_genomes\treference_genomes\tdiscarded_genomes\n')
        fl.write(f'{args.species}\t{selection_no_ref}\t{reference_genome_count}\t{discarded_genome_count}\n')


if __name__ == "__main__":
    main()
