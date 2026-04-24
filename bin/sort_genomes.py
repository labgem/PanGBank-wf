#!/usr/bin/env python

import argparse
import logging
import sys
from pathlib import Path
import numpy as np
import pandas as pd


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--genome_name_to_path",
        type=Path,
        required=True,
        help="TSV with two columns: genome name and genome path. No header.",
    )

    parser.add_argument(
        "--genome_stats",
        type=Path,
        required=True,
        help="Path to the input TSV file containing genome statistics, which will be used for sorting.",
    )

    parser.add_argument(
        "--genome_metadata",
        type=Path,
        required=False,
        help="Path to the input TSV file containing genome metadata.",
    )

    parser.add_argument(
        "--completeness_sorting",
        action="store_true",
        default=False,
        help=(
            "Enable sorting based on metadata completeness metrics in addition to assembly statistics. "
            "Requires --genome_metadata and the columns specified by --completeness_columns to be present."
        ),
    )

    parser.add_argument(
        "--sorted_genome_list",
        type=Path,
        default="sorted_genomes.txt",
        help="Path to the output file where the sorted list of genome identifiers will be written.",
    )

    parser.add_argument(
        "--sorted_genome_stats",
        type=Path,
        default="sorted_genomes.tsv",
        help="Path to the output TSV file containing the sorted genome statistics.",
    )

    parser.add_argument(
        "-l",
        "--log-level",
        help="Desired log level (default: INFO).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="INFO",
    )

    return parser.parse_args(argv)


STATS_SORT_COLUMNS = ["L90", "L75", "L50", "auN"]
STATS_ASCENDING_COLUMNS = ["L50", "L75", "L90"]

COMPLETENESS_COLUMNS = ["checkm2_completeness", "checkm_completeness"]

# Extra columns from metadata kept in the output (in addition to all stats columns).
# Missing ones trigger a warning but not an error.
OUTPUT_EXTRA_COLUMNS = [
    "checkm2_completeness",
    "checkm2_contamination",
    "checkm_completeness",
    "checkm_contamination",
    "max_completeness",
]

def select_output_columns(df, stats_columns):
    """Return df restricted to stats columns + genome_name + available OUTPUT_EXTRA_COLUMNS.

    Columns listed in OUTPUT_EXTRA_COLUMNS that are absent from df trigger a warning.
    """
    extra = []
    for col in OUTPUT_EXTRA_COLUMNS:
        if col in df.columns:
            extra.append(col)
        else:
            logging.warning(f"Output column '{col}' is not present in the data and will be skipped.")
    return df[[c for c in stats_columns + ["genome_name"] + extra if c in df.columns]]


def check_no_nulls(df, columns, context):
    """Raise ValueError if any of the given columns contain NaN values."""
    for col in columns:
        null_genomes = df.loc[df[col].isna(), "genome_name"].tolist()
        if null_genomes:
            raise ValueError(
                f"Column '{col}' ({context}) has missing values for "
                f"{len(null_genomes)} genome(s): {', '.join(str(g) for g in null_genomes)}."
            )

def sort_by_stats(df):
    """Sort genomes based on assembly statistics only.

    L50/L75/L90 sorted ascending (lower index = better);
    auN sorted descending (larger value = better).
    """
    ascending_sorts = [col in STATS_ASCENDING_COLUMNS for col in STATS_SORT_COLUMNS]
    logging.info(
        f"Sorting by assembly statistics: {', '.join(STATS_SORT_COLUMNS)} "
        f"(ascending: {[c for c, a in zip(STATS_SORT_COLUMNS, ascending_sorts) if a]}, "
        f"descending: {[c for c, a in zip(STATS_SORT_COLUMNS, ascending_sorts) if not a]})"
    )
    return df.sort_values(by=STATS_SORT_COLUMNS, ascending=ascending_sorts)


def sort_by_completeness(df):
    """Sort genomes by floored max completeness across CheckM/CheckM2, then by assembly statistics.

    Step 1 — derive max_completeness (integer, descending):
        Take the highest completeness value available across checkm_completeness and
        checkm2_completeness for each genome, then floor it to the nearest integer.
        This groups genomes into 1 % completeness bins so that small floating-point
        differences between tools do not dominate the ranking.

    Step 2 — break ties with assembly statistics (in order):
        L90  ascending  — fewer contigs needed to cover 90 % of assembly = better
        L75  ascending  — fewer contigs needed to cover 75 % of assembly = better
        L50  ascending  — fewer contigs needed to cover 50 % of assembly = better
        auN  descending — larger area-under-Nx-curve = better contiguity
    """
    logging.info("Computing max_completeness: floor(max(checkm2_completeness, checkm_completeness))")
    df = df.copy()
    df["max_completeness"] = (
        np.floor(df[COMPLETENESS_COLUMNS].max(axis=1, skipna=True))
        .astype("Int64")
    )

    sort_columns = ["max_completeness"] + STATS_SORT_COLUMNS
    ascending_sorts = [col in STATS_ASCENDING_COLUMNS for col in sort_columns]

    logging.info(
        "Sort order: "
        + ", ".join(
            f"{col} ({'asc' if asc else 'desc'})"
            for col, asc in zip(sort_columns, ascending_sorts)
        )
    )

    return df.sort_values(by=sort_columns, ascending=ascending_sorts)


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)

    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")

    # --- Input validation ---
    if not args.genome_stats.is_file():
        raise FileNotFoundError(f"Genome stat file {args.genome_stats} was not found!")

    if not args.genome_name_to_path.is_file():
        raise FileNotFoundError(f"genome_name_to_path file {args.genome_name_to_path} was not found!")

    if not args.sorted_genome_list.parent.exists():
        raise FileNotFoundError(
            f"Output directory for sorted_genome_list {args.sorted_genome_list.parent} was not found!"
        )

    if args.sorted_genome_stats and not args.sorted_genome_stats.parent.exists():
        raise FileNotFoundError(
            f"Output directory for sorted_genome_stats {args.sorted_genome_stats.parent} was not found!"
        )

    if args.completeness_sorting and not args.genome_metadata:
        raise ValueError("--completeness_sorting requires --genome_metadata to be provided.")

    if args.genome_metadata and not args.genome_metadata.is_file():
        raise FileNotFoundError(f"genome_metadata file {args.genome_metadata} was not found!")

    # --- Load data ---
    logging.info(f"Loading genome statistics from {args.genome_stats}")
    df = pd.read_csv(args.genome_stats, sep="\t")
    stats_columns = df.columns.tolist()  # record original stats columns before any merge

    unfound_stat_columns = [col for col in STATS_SORT_COLUMNS if col not in df.columns]
    if unfound_stat_columns:
        raise ValueError(
            f"The following required stat columns were not found in {args.genome_stats}: "
            f"{', '.join(unfound_stat_columns)}. "
            f"Available columns: {', '.join(df.columns)}."
        )

    logging.info(f"Loading genome name-to-path mapping from {args.genome_name_to_path}")
    name_to_path_df = pd.read_csv(
        args.genome_name_to_path, sep="\t", header=None, names=["genome_name", "genome_path"]
    )
    df = df.merge(name_to_path_df, left_on="File", right_on="genome_path", how="left")
    df = df.drop(columns=["genome_path"])

    if args.genome_metadata:
        logging.info(f"Loading genome metadata from {args.genome_metadata}")
        metadata_df = pd.read_csv(args.genome_metadata, sep="\t")

        if args.completeness_sorting:
            unfound_meta_columns = [col for col in COMPLETENESS_COLUMNS if col not in metadata_df.columns]
            if unfound_meta_columns:
                raise ValueError(
                    f"The following required completeness columns were not found in {args.genome_metadata}: "
                    f"{', '.join(unfound_meta_columns)}. "
                    f"Available columns: {', '.join(metadata_df.columns)}."
                )

        df = df.merge(metadata_df, left_on="genome_name", right_on="genomes", how="left")
        df = df.drop(columns=["genomes"])

    # --- Sort ---
    if args.completeness_sorting:
        check_no_nulls(df, COMPLETENESS_COLUMNS, "completeness sorting")
        check_no_nulls(df, STATS_SORT_COLUMNS, "assembly statistics")
        sorted_df = sort_by_completeness(df)
    else:
        logging.info("Completeness sorting disabled — sorting on assembly statistics only.")
        check_no_nulls(df, STATS_SORT_COLUMNS, "assembly statistics")
        sorted_df = sort_by_stats(df)

    # --- Output ---
    if args.sorted_genome_stats:
        logging.info(f"Writing sorted genome stats to {args.sorted_genome_stats}")
        select_output_columns(sorted_df, stats_columns).to_csv(args.sorted_genome_stats, sep="\t", index=False)

    logging.info(f"Writing sorted genome list to {args.sorted_genome_list}")
    sorted_df["File"].to_csv(args.sorted_genome_list, sep="\t", index=False, header=False)


if __name__ == "__main__":
    sys.exit(main())
