#!/usr/bin/env python

import argparse
import csv
import logging
import sys
from pathlib import Path
import yaml


def get_info_from_yaml_pangenome_info(yaml_info):

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


def get_info_from_yaml_genome_stat_summary(yaml_info):

    name = yaml_info.stem

    with open(yaml_info, "r") as fh:
        pangenome_info = yaml.safe_load(fh)

    metric_of_interest = [
        "Persistent_families",
        "Soft_core_families",
        "Exact_core_families",
        "Shell_families",
        "Cloud_families",
    ]

    useful_info = {"Name": name}
    for metric in metric_of_interest:
        for stat in ["mean"]:
            useful_info[f"{metric}_{stat}".capitalize()] = pangenome_info[metric][stat]

    return useful_info


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

    parser.add_argument(
        "--yaml_info_dir",
        help="Directory where yaml info files are stored",
        required=True,
        type=Path,
    )
    parser.add_argument(
        "--yaml_genome_stats_dir",
        help="Directory where summarize genome stat files are stored",
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

    return args


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)

    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")

    name_to_info = {}

    for i, yaml_info in enumerate(args.yaml_info_dir.iterdir()):
        logging.info(f"{i}: {yaml_info}")

        info = get_info_from_yaml_pangenome_info(yaml_info)
        name_to_info[info["Name"]] = info

    if args.yaml_genome_stats_dir:

        for i, yaml_genome_stat_summary in enumerate(
            args.yaml_genome_stats_dir.iterdir()
        ):
            logging.info(f"{i}: {yaml_info}")

            info = get_info_from_yaml_genome_stat_summary(yaml_genome_stat_summary)
            name_to_info[info["Name"]].update(info)

    info_to_write = list(name_to_info.values())

    write_tsv_from_list_of_dict(args.output, info_to_write)


if __name__ == "__main__":
    sys.exit(main())
