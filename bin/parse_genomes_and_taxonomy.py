#!/usr/bin/env python


"""Provide a command line tool to parse genome paths and taxonomy."""


import argparse
import csv
import logging
import sys
from collections import defaultdict
from pathlib import Path
import gzip
import re


# TODO check circular contig when input genomes are in fasta and add them in ppanggolin input files


def parse_genome_files(genomes_paths_file):
    """ """
    files_not_found = []

    acc_to_genome_file = {}
    accession_count = 0

    with open(genomes_paths_file) as fl:
        for i, line in enumerate(fl):
            if not line:
                continue

            genome_file = Path(line.strip().split()[1])
            name = line.strip().split()[0]

            if not genome_file.is_file():
                files_not_found.append((i + 1, name, genome_file))

            accession_count += 1
            acc_to_genome_file[name] = genome_file

    assert accession_count == len(acc_to_genome_file), "Some genome names are duplicated in the genome path file."

    if files_not_found:
        for line, acc, genome_file in files_not_found:
            logging.error(f"Genome file '{genome_file}' at line {line} of the genome path file '{genomes_paths_file}' was not found!")
        sys.exit(2)

    return acc_to_genome_file


def parse_taxonomy_file(taxonomy_file):
    """ """

    proper_open = gzip.open if taxonomy_file.suffix == ".gz" else open

    sptax_to_accessions = defaultdict(set)
    # acc_to_sptax = {}

    with proper_open(taxonomy_file, "rt") as fl:
        for line in fl:
            accession, taxonomy = line.strip().split("\t")

            if accession.startswith("RS_") or accession.startswith("GB_"):
                accession = "_".join(accession.split("_")[1:])

            sptax_to_accessions[taxonomy].add(accession)
            # acc_to_sptax[accession] = taxonomy

    return sptax_to_accessions


def associate_genomes_and_taxonomy(acc_to_genome_file, sptax_to_accessions, min_genome_count):
    """ """
    remove_assembly_version = lambda accessions: {acc.split(".")[0] for acc in accessions}

    all_input_accessions = remove_assembly_version(acc_to_genome_file)

    acc_to_acc_with_version = {acc.split(".")[0]: acc for acc in acc_to_genome_file}

    sptax_to_input_accs = {}

    species_infos = []
    for sptax, sp_accessions in sptax_to_accessions.items():
        input_sp_accessions = all_input_accessions & remove_assembly_version(sp_accessions)

        sptax_to_input_accs[sptax] = {acc_to_acc_with_version[acc] for acc in input_sp_accessions}

        species = sptax.split(";")[-1]

        sp_genbank_accesions_count = len([acc for acc in sp_accessions if acc.startswith("GCA_")])
        sp_refseq_accesions_count = len([acc for acc in sp_accessions if acc.startswith("GCF_")])

        build_pangenome = True if len(input_sp_accessions) >= min_genome_count else False

        sp_info = {
            "species": species,
            "sp_genome_in_taxonomy": len(sp_accessions),
            "sp_genome_from_refseq": sp_refseq_accesions_count,
            "sp_genome_from_genbank": sp_genbank_accesions_count,
            "input_sp_genome": len(input_sp_accessions),
            "build_pangenome": build_pangenome,
            "Taxonomy": sptax,
        }

        species_infos.append(sp_info)

    return sptax_to_input_accs, species_infos


def write_species_summary(species_summary_file, species_infos):
    """ """
    # Extract column names from the keys of the first dictionary
    fieldnames = species_infos[0].keys()

    with open(species_summary_file, "w", newline="", encoding="utf-8") as tsvfile:
        writer = csv.DictWriter(tsvfile, fieldnames=fieldnames, delimiter="\t")

        writer.writeheader()

        writer.writerows(species_infos)


def write_ppanggolin_input_files(outdir, sptax_to_input_accs_filtered, acc_to_genome_file):
    """ """
    for sptax, input_accs in sptax_to_input_accs_filtered.items():
        species = sptax.split(";")[-1].replace(" ", "_")

        with open(outdir / f"{species}.tsv", "w") as flout:
            flout.write("\n".join(f"{acc}\t{acc_to_genome_file[acc]}" for acc in input_accs) + "\n")


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Parse genome paths and taxonomy.",
        epilog="Example: python parse_genomes_and_taxonomy.py --genomes <genomes_file> --taxonomy <taxonomy_file>",
    )
    parser.add_argument(
        "--genomes",
        type=Path,
        required=True,
        help="Path to a TSV file containing input genome paths."
        "The file is expected to have two columns: the first column containing genome accessions, "
        "and the second column containing the path to genome file.",
    )
    parser.add_argument(
        "--taxonomy",
        type=Path,
        required=True,
        help="Path to the taxonomy file in TSV format corresponding to the input genomes. "
        "The file is expected to have two columns: the first column containing genome accessions, "
        "and the second column containing taxonomy information separated by semicolons. "
        "The last taxon of the taxonomy is expected to represent the species.",
    )

    parser.add_argument(
        "--min_genomes",
        help="Minimum number of genomes required to build a pangenome.",
        default=15,
        type=int,
    )

    parser.add_argument(
        "-o",
        "--ppanggolin_files_outdir",
        help="Directory where ppanggolin input files will be created.",
        default="ppanggolin_input_files",
        type=Path,
    )

    parser.add_argument(
        "--species_summary_file",
        help="Output file containing species summary information, including the number of genomes and whether they meet the pangenome threshold.",
        default="species_summary.tsv",
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

    if not args.genomes.is_file():
        logging.error(f"Genome paths file {args.genomes} was not found!")
        sys.exit(2)

    if not args.taxonomy.is_file():
        logging.error(f"Taxonomy file {args.taxonomy} was not found!")
        sys.exit(2)

    if not args.species_summary_file.parent.exists():
        logging.error(f"The directory of species_summary_outfile {args.species_summary_file.parent} was not found!")
        sys.exit(2)

    args.ppanggolin_files_outdir.mkdir(parents=True, exist_ok=True)

    acc_to_genome_file = parse_genome_files(args.genomes)

    sptax_to_accessions = parse_taxonomy_file(args.taxonomy)

    sptax_to_input_accs, species_infos = associate_genomes_and_taxonomy(acc_to_genome_file, sptax_to_accessions, args.min_genomes)

    logging.info(f"Writing species summary in {args.species_summary_file}")
    write_species_summary(args.species_summary_file, species_infos)

    sptax_to_input_accs_filtered = {sptax: accs for sptax, accs in sptax_to_input_accs.items() if len(accs) >= args.min_genomes}

    logging.info(f"{len(sptax_to_input_accs_filtered)} species have enough genomes to build a pangenome.")

    logging.info(f"Writing ppanggolin input files in {args.ppanggolin_files_outdir}")
    write_ppanggolin_input_files(args.ppanggolin_files_outdir, sptax_to_input_accs_filtered, acc_to_genome_file)


if __name__ == "__main__":
    sys.exit(main())
