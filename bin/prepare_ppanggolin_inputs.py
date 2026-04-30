#!/usr/bin/env python


"""Provide a command line tool to parse genome paths and taxonomy."""


import argparse
import logging
import sys
from pathlib import Path
import pandas as pd

# TODO check circular contig when input genomes are in fasta and add them in ppanggolin input files


def parse_genome_files(genomes_paths_file, check_files_existence=True) -> pd.DataFrame:
    """Return a DataFrame with columns ['name', 'path']."""
    logging.info(f"Parsing genome file {genomes_paths_file}")
    df = pd.read_csv(genomes_paths_file, sep="\t", names=["name", "path"], dtype=str)

    dup_names = df.loc[df.duplicated("name", keep=False), "name"].unique()
    if len(dup_names):
        raise ValueError(
            f"Duplicate genome names in {genomes_paths_file}: {', '.join(dup_names)}"
        )

    if check_files_existence:
        missing = df[~df["path"].map(lambda p: Path(p).is_file())]
        if not missing.empty:
            for row in missing.itertuples():
                logging.error(
                    f"Genome file '{row.path}' at line {row.Index + 1} of '{genomes_paths_file}' was not found!"
                )
            sys.exit(2)

    return df


def parse_taxonomy_file(taxonomy_file) -> pd.DataFrame:
    """Return a DataFrame with columns ['accession', 'taxonomy']."""
    logging.info(f"Parsing taxonomy file {taxonomy_file}")
    return pd.read_csv(
        taxonomy_file,
        sep="\t",
        names=["accession", "taxonomy"],
        dtype=str,
    )


def associate_genomes_and_taxonomy(
    genomes_df: pd.DataFrame, taxonomy_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Match input genomes to their taxonomy.

    Returns:
        matched_df: one row per genome with columns ['name', 'path', 'taxonomy']
    """
    # Strip assembly version suffix so GCA_X.1 and GCA_X.2 match the same genome
    tax = taxonomy_df.copy()
    tax["acc_no_version"] = tax["accession"].str.rsplit(".", n=1).str[0]

    genomes = genomes_df.copy()
    genomes["acc_no_version"] = genomes["name"].str.rsplit(".", n=1).str[0]

    merged = tax.merge(genomes, on="acc_no_version", how="inner")

    missing = genomes.loc[~genomes["name"].isin(merged["name"]), "name"]
    if not missing.empty:
        raise ValueError(
            f"{len(missing)} genome(s) from the genome file have no matching taxonomy:\n"
            + "\n".join(f"  - {name}" for name in missing)
        )

    dup_names = merged.loc[merged.duplicated("name", keep=False), "name"].unique()
    if len(dup_names):
        raise ValueError(
            f"{len(dup_names)} genome(s) matched multiple taxonomy entries:\n"
            + "\n".join(f"  - {name}" for name in dup_names)
        )

    matched_df = merged[["name", "path", "taxonomy"]].copy()

    # Ensure no species maps to more than one taxonomy string
    species_col = matched_df["taxonomy"].str.split(";").str[-1]
    multi = (
        matched_df.groupby(species_col)["taxonomy"]
        .apply(set)
        .loc[lambda s: s.map(len) > 1]
    )
    if not multi.empty:
        for sp, taxs in multi.items():
            logging.warning(f"Species {sp} has multiple taxonomies: {taxs}")
        raise ValueError(f"{len(multi)} species have multiple taxonomies.")

    return matched_df


def filter_input_genome(genome_df: pd.DataFrame, genomes: set[str]) -> pd.DataFrame:
    filtered = genome_df[genome_df["genomes"].isin(genomes)]
    logging.info(
        f"Restricting taxonomy to input genomes: {len(filtered)}/{len(genome_df)} entries retained "
        f"({len(genome_df) - len(filtered)} not found in taxonomy)"
    )
    return filtered


def write_species_summary(
    species_summary_file,
    matched_df: pd.DataFrame,
    taxonomy_df: pd.DataFrame,
    min_genome_count: int,
):
    """Build a per-taxonomy summary and write it to species_summary_file."""
    logging.info(f"Writing species summary in {species_summary_file}")

    sp_total = (
        taxonomy_df.groupby("taxonomy")["accession"]
        .nunique()
        .rename("sp_genome_in_taxonomy")
    )
    sp_input = (
        matched_df.groupby("taxonomy")["name"].nunique().rename("input_sp_genome")
    )

    summary_df = sp_total.reset_index().merge(
        sp_input.reset_index(), on="taxonomy", how="left"
    )
    summary_df["input_sp_genome"] = summary_df["input_sp_genome"].fillna(0).astype(int)
    summary_df["species"] = summary_df["taxonomy"].str.split(";").str[-1]
    summary_df["build_pangenome"] = summary_df["input_sp_genome"] >= min_genome_count
    summary_df = summary_df[
        [
            "species",
            "sp_genome_in_taxonomy",
            "input_sp_genome",
            "build_pangenome",
            "taxonomy",
        ]
    ]

    summary_df.to_csv(species_summary_file, sep="\t", index=False)


def filter_genomes_for_pangenome(
    matched_df: pd.DataFrame, min_genome_count: int, uninformative_species: set[str]
) -> pd.DataFrame:
    """Return genomes belonging to species with enough genomes and an informative taxonomy."""
    sp_counts = matched_df.groupby("taxonomy")["name"].transform("count")
    filtered_df = matched_df[
        (sp_counts >= min_genome_count)
        & ~matched_df["taxonomy"].str.split(";").str[-1].isin(uninformative_species)
    ].copy()

    logging.info(
        f"{filtered_df['taxonomy'].nunique()} species have enough genomes to build a pangenome."
    )
    return filtered_df


def write_ppanggolin_input_files(outdir, filtered_df: pd.DataFrame):
    logging.info(f"Writing ppanggolin input files in {outdir}")
    for taxonomy, group in filtered_df.groupby("taxonomy"):
        species = taxonomy.split(";")[-1].replace(" ", "_")
        sp_outdir = outdir / species
        sp_outdir.mkdir(parents=True, exist_ok=True)
        group[["name", "path"]].sort_values("name").to_csv(
            sp_outdir / "input_genomes.tsv.gz",
            sep="\t",
            index=False,
            header=False,
            compression="gzip",
        )


def write_metadata_by_species(outdir, genome_metadata_df: pd.DataFrame):
    logging.info(f"Splitting and writing genome metadata by species in {outdir}")
    for sptax, group in genome_metadata_df.groupby("_sptax"):

        species = sptax.split(";")[-1].replace(" ", "_")
        sp_outdir = outdir / species
        sp_outdir.mkdir(parents=True, exist_ok=True)

        group.drop(columns=["_sptax"]).to_csv(
            sp_outdir / "genomes_metadata.tsv.gz",
            sep="\t",
            index=False,
            compression="gzip",
        )


def parse_metadata_file(
    genome_metadata_file: Path, taxonomy_map: pd.Series
) -> pd.DataFrame:
    df = pd.read_csv(genome_metadata_file, sep="\t", dtype=str, keep_default_na=False)

    if "genomes" not in df.columns:
        raise KeyError(
            f"genomes column not found in metadata file {genome_metadata_file}"
        )

    df = df[df["genomes"].isin(taxonomy_map.index)].copy()
    df["_sptax"] = df["genomes"].map(taxonomy_map)

    return df


def parse_translation_table_file(
    translation_table_file: Path, taxonomy_map: pd.Series
) -> pd.DataFrame:
    """
    Parse translation table file and associate translation tables with species taxonomy.

    Args:
        translation_table_file: Path to TSV file containing genome accessions and translation tables
        taxonomy_map: Series indexed by genome name, mapping to taxonomy string

    Returns:
        DataFrame with columns ['name', 'translation_table', 'taxonomy']
    """
    logging.info(f"Parsing translation table file {translation_table_file}")
    df = pd.read_csv(
        translation_table_file,
        sep="\t",
        header=None,
        names=["name", "translation_table"],
        dtype=str,
        keep_default_na=False,
    )

    invalid_mask = ~df["translation_table"].str.isdigit()
    for name in df.loc[invalid_mask, "name"]:
        logging.debug(f"Invalid translation table value for genome {name}. Skipping.")
    df = df[~invalid_mask]

    df = df[df["name"].isin(taxonomy_map.index)].copy()
    df["taxonomy"] = df["name"].map(taxonomy_map)

    logging.info(
        f"Parsed translation table for {len(df)} genomes across {df['taxonomy'].nunique()} species."
    )
    return df


def write_translation_table_by_species(outdir, translation_df: pd.DataFrame):
    logging.info(f"Splitting and writing translation table by species in {outdir}")
    rows = []
    for taxonomy, group in translation_df.groupby("taxonomy"):
        species = taxonomy.split(";")[-1].replace(" ", "_")
        counts = group["translation_table"].value_counts()
        if len(counts) > 1:
            logging.warning(
                f"Species {species} has multiple translation tables: {counts.to_dict()}. The most common one will be used."
            )
        rows.append(f"{species}\t{counts.index[0]}\n")
    (outdir / "species_to_translation_tables.tsv").write_text("".join(rows))


def existing_file(value: str) -> Path:
    """Argparse type: return Path only if the file exists."""
    path = Path(value)
    if not path.is_file():
        raise argparse.ArgumentTypeError(f"file not found: {value}")
    return path


def existing_parent(value: str) -> Path:
    """Argparse type: return Path only if the parent directory exists."""
    path = Path(value)
    if not path.parent.exists():
        raise argparse.ArgumentTypeError(
            f"parent directory does not exist: {path.parent}"
        )
    return path


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Parse genome paths and taxonomy.",
        epilog="Example: python parse_genomes_and_taxonomy.py --genomes <genomes_file> --taxonomy <taxonomy_file>",
    )
    parser.add_argument(
        "--genomes",
        type=existing_file,
        required=True,
        help="Path to a TSV file containing input genome paths."
        "The file is expected to have two columns: the first column containing genome accessions, "
        "and the second column containing the path to genome file.",
    )
    parser.add_argument(
        "--taxonomy",
        type=existing_file,
        required=True,
        help="Path to the taxonomy file in TSV format corresponding to the input genomes. "
        "The file is expected to have two columns: the first column containing genome accessions, "
        "and the second column containing taxonomy information separated by semicolons. "
        "The last taxon of the taxonomy is expected to represent the species.",
    )

    parser.add_argument(
        "--species_to_merge",
        type=existing_file,
        required=False,
        help="Path to the file containing species to merge information in TSV format. "
        "The file is expected to have two columns: the first column containing species identifiers, "
        "and the second column containing the species to merge.",
    )

    parser.add_argument(
        "--genome_metadata",
        type=existing_file,
        required=False,
        help="Path to a metadata file in TSV format corresponding to the input genomes. "
        "Expected to have a column genomes containing genome accessions.",
    )

    parser.add_argument(
        "--genome_translation_table",
        type=existing_file,
        required=False,
        help="Path with genome accessions and their corresponding translation table in TSV format. ",
    )

    parser.add_argument(
        "--min_genomes",
        help="Minimum number of genomes required to build a pangenome.",
        default=15,
        type=int,
    )

    parser.add_argument(
        "-o",
        "--outdir",
        help="Directory where species dir will be create to store ppanggolin input files and genome metadata by species if provided.",
        default="ppanggolin_input_files",
        type=Path,
    )

    parser.add_argument(
        "--species_summary_file",
        help="Output file containing species summary information, including the number of genomes and whether they meet the pangenome threshold.",
        default="species_summary.tsv",
        type=existing_parent,
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

    uninformative_species = {"s__"}  # useful for mgnify genomes taxonomy

    args.outdir.mkdir(parents=True, exist_ok=True)

    genome_path_df = parse_genome_files(args.genomes, check_files_existence=False)
    genome_taxonomy_df = parse_taxonomy_file(args.taxonomy)

    matched_df = associate_genomes_and_taxonomy(genome_path_df, genome_taxonomy_df)

    write_species_summary(
        args.species_summary_file, matched_df, genome_taxonomy_df, args.min_genomes
    )

    # Filter to species suitable for pangenome building
    filtered_df = filter_genomes_for_pangenome(
        matched_df, args.min_genomes, uninformative_species
    )

    write_ppanggolin_input_files(args.outdir, filtered_df)

    taxonomy_map = filtered_df.set_index("name")["taxonomy"]

    if args.genome_metadata:
        metadata_df = parse_metadata_file(args.genome_metadata, taxonomy_map)
        write_metadata_by_species(args.outdir, metadata_df)

    if args.genome_translation_table:
        translation_df = parse_translation_table_file(
            args.genome_translation_table, taxonomy_map
        )
        write_translation_table_by_species(args.outdir, translation_df)

if __name__ == "__main__":
    sys.exit(main())
