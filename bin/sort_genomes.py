#!/usr/bin/env python

import argparse
import logging
import sys
from pathlib import Path
import pandas as pd


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
    )
    parser.add_argument(
        "--genome_stats",
        type=Path,
        required=True,
        help="Path to the input TSV file containing genome statistics, which will be used for sorting.",
    )

    parser.add_argument(
        "--sorted_genome_list",
        type=Path,
        default="sorted_genomes.txt",
        help=(
            "Path to the output file where the sorted list of genome identifiers will be written."
        ),
    )

    parser.add_argument(
        "--sorted_genome_stats",
        type=Path,
        default="sorted_genomes.tsv",
        help=(
            "Path to the output TSV file containing the sorted genome statistics. "
            "This file will be based on the input genome statistics. Default is to not write this file."
        ),
    )


    parser.add_argument(
        "--sort_by",
        type=str,
        nargs='+',
        default=['L90', 'L75', 'L50', 'auN'],
        help=(
            "Specify the columns to sort genomes by. Sorting is ascending for 'genome index' metrics "
            "(L50, L75, L90) and descending for all other columns that represent contig sizes "
            "(e.g., N50, auN, Total bp)."
        ),
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

    if not args.genome_stats.is_file():
        logging.error(f"Genome stat file {args.genome_stats} was not found!")
        sys.exit(2)

    if not args.sorted_genome_list.parent.exists():
        logging.error(f"The directory of sorted_genomes output file {args.sorted_genomes.parent} was not found!")
        sys.exit(2)

    if args.sorted_genome_stats and not args.sorted_genome_stats.parent.exists():
        logging.error(f"The directory of sorted_genomes output file {args.sorted_genomes.parent} was not found!")
        sys.exit(2)

    logging.info(f"Parsing genome file {args.genome_stats}")

    df = pd.read_csv(args.genome_stats, sep="\t")

    columns_to_sort = args.sort_by  # Renamed for better clarity
    # Check if the columns exist
    unfound_columns = [column for column in columns_to_sort if column not in df.columns]

    if unfound_columns:
        raise ValueError(
            f"The following columns were not found in the genome stats file ({args.genome_stats}): "
            f"{', '.join(unfound_columns)}.\n"
            f"Available columns are: {', '.join(df.columns)}."
        )

    logging.info(f"Sorting genome based on {columns_to_sort}")

    ascending_columns = ['L50', 'L75','L90']
    ascending_sorts = [column in ascending_columns for column in columns_to_sort]
    # Sort by the specified column
    sorted_df = df.sort_values(by=columns_to_sort, ascending=ascending_sorts)

    # Output the result
    if args.sorted_genome_stats:
        logging.info(f"Writing sorted genome stats in {args.sorted_genome_stats}")
        sorted_df.to_csv(args.sorted_genome_stats, sep="\t", index=False)


    logging.info(f"Writing sorted genome list in {args.sorted_genome_list}")
    sorted_df['File'].to_csv(args.sorted_genome_list, sep="\t", index=False, header=False)


if __name__ == "__main__":
    sys.exit(main())
