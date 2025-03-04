#!/usr/bin/env python

import argparse
import csv
import logging
import sys
from pathlib import Path
import yaml
import pandas as pd


def get_info_from_yaml(yaml_info):

    name = yaml_info.stem

    with open(yaml_info, "r") as fh:
        pangenome_info = yaml.safe_load(fh)["Content"]

    useful_info = {
        "Name": name,
        "Genomes": pangenome_info["Genomes"],
        "Genes": pangenome_info["Genes"],
        "Families": pangenome_info["Families"],
        "Persistent": pangenome_info["Persistent"]["Family_count"],
        "Shell": pangenome_info["Shell"]["Family_count"],
        "Cloud": pangenome_info["Cloud"]["Family_count"],
        "RGPs": pangenome_info["RGP"],
        "Spots": pangenome_info["Spots"],
        "Modules": pangenome_info["Modules"]["Number_of_modules"],
        "Partitions": pangenome_info["Number_of_partitions"],
    }
    return useful_info


def summarize_genome_stat(genome_stat_file):
    """ """

    df = pd.read_csv(
        genome_stat_file,
        sep="\t",
        comment="#",
    )

    df["Persistent_families_fraction"] = df["Persistent_families"] / df["Families"]
    df["Soft_core_families_fraction"] = df["Soft_core_families"] / df["Families"]
    df["Exact_core_families_fraction"] = df["Exact_core_families"] / df["Families"]

    df["Shell_families_fraction"] = df["Shell_families"] / df["Families"]
    df["Shell_families_fraction"] = df["Cloud_families"] / df["Families"]

    df["Variable_families"] = df["Shell_families"] + df["Cloud_families"]

    df["Variable_families_fraction"] = df["Variable_families"] / df["Families"]

    columns_to_process = [
        "Persistent_families_fraction",
        "Soft_core_families_fraction",
        "Exact_core_families_fraction",
        # "Shell_families_fraction",
        # "Variable_families_fraction",
        "Fragmentation",
        "Completeness",
        "Contamination",
        # "Genes",
        # "Contigs",
    ]

    species_stats = {}

    # Calculate stats for each column
    # operations = ["min", "max", "mean", "median", "std"]
    operations = ["median"]
    for column in columns_to_process:
        if column in df.columns:
            stats = df[column].agg(operations).to_dict()

            stats["Q1"] = df[column].quantile(0.25)
            stats["Q3"] = df[column].quantile(0.75)

            # Compute IQR
            stats["IQR"] = stats["Q3"] - stats["Q1"]

            for stat_name, value in stats.items():
                species_stats[f"{stat_name}_{column}"] = value

    return species_stats


def write_tsv_from_list_of_dict(species_summary_file, species_infos):
    """ """

    # Extract column names from the keys of the first dictionary
    fieldnames = species_infos[0].keys()

    with open(species_summary_file, "w", newline="", encoding="utf-8") as tsvfile:
        writer = csv.DictWriter(tsvfile, fieldnames=fieldnames, delimiter="\t")

        writer.writeheader()

        writer.writerows(species_infos)

def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Parse genome paths and taxonomy.",
        epilog="Example: python gather_pangenome_infos.py --yaml_dir <yaml_info_dir>",
    )

    group1 = parser.add_argument_group("Single file arguments")
    group1.add_argument(
        "--yaml_info",
        help="Directory where yaml info files are stored",
        required=False,
        type=Path,
    )
    group1.add_argument(
        "--genome_stat",
        help="Directory where genome stat files are stored",
        required=False,
        type=Path,
    )

    group2 = parser.add_argument_group("Directory arguments")
    group2.add_argument(
        "--yaml_dir",
        help="Directory where yaml info files are stored",
        required=False,
        type=Path,
    )
    group2.add_argument(
        "--genome_stat_dir",
        help="Directory where genome stat files are stored",
        required=False,
        type=Path,
    )

    parser.add_argument(
        "-o",
        "--output",
        help="Output file containing pangenome infos.",
        default="pangenome_summary.tsv",
        type=Path,
    )

    parser.add_argument(
        "-l",
        "--log-level",
        help="Desired log level (default: WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="INFO",
    )

    args = parser.parse_args(argv)

    # Ensure mutual exclusivity between (yaml_info + genome_stat) and (yaml_dir + genome_stat_dir)
    single_file_args = args.yaml_info is not None or args.genome_stat is not None
    dir_args = args.yaml_dir is not None or args.genome_stat_dir is not None

    if single_file_args and dir_args:
        parser.error(
            "Cannot mix --yaml_info/--genome_stat with --yaml_dir/--genome_stat_dir. Choose one set."
        )

    if not single_file_args and not dir_args:
        parser.error(
            "You must provide either (--yaml_info and --genome_stat) or (--yaml_dir and --genome_stat_dir)."
        )

    return args


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)

    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")

    name_to_info = {}

    if args.yaml_dir:
        for i, yaml_info in enumerate(args.yaml_dir.iterdir()):
            logging.info(f"{i}: {yaml_info}")

            info = get_info_from_yaml(yaml_info)
            name_to_info[info["Name"]] = info

        if args.genome_stat_dir:
            for i, genome_stat_dir in enumerate(args.genome_stat_dir.iterdir()):
                if not genome_stat_dir.is_dir():
                    continue
                genome_stat_file = genome_stat_dir / "genomes_statistics.tsv.gz"
                name = genome_stat_dir.name
                logging.info(f"{i}: {genome_stat_file}")

                genome_stat_summary = summarize_genome_stat(genome_stat_file)

                info = name_to_info[name]

                info.update(genome_stat_summary)

        info_to_write = list(name_to_info.values())

    if args.yaml_info:
        info = get_info_from_yaml(args.yaml_info)
        if args.genome_stat:
            genome_stat_summary = summarize_genome_stat(args.genome_stat)
            info.update(genome_stat_summary)

        info_to_write = [info]

    write_tsv_from_list_of_dict(args.output, info_to_write)


if __name__ == "__main__":
    sys.exit(main())
