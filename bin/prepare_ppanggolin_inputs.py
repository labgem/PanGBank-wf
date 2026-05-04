#!/usr/bin/env python


"""Provide a command line tool to parse genome paths and taxonomy."""


import argparse
import logging
import re
import sys
from pathlib import Path
import pandas as pd


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
                    f"Genome file '{row.path}' at line {row.Index + 1} of '{genomes_paths_file}' was not found!"  # type: ignore
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


def parse_species_to_merge_file(species_to_merge_file, min_genomes) -> pd.DataFrame:
    """Return a DataFrame with columns ['meta_sp_dir_name', 'species_to_merge']."""
    logging.info(f"Parsing species to merge file {species_to_merge_file}")
    df = pd.read_csv(
        species_to_merge_file,
        sep="\t",
        names=["meta_sp_dir_name", "species_to_merge", "genome_count"],
        dtype=str,
    )
    df = df.loc[df["genome_count"].astype(int) >= min_genomes]
    return df


def associate_genomes_and_taxonomy(
    genomes_df: pd.DataFrame, taxonomy_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Match input genomes to their taxonomy.

    Returns:
        matched_df: one row per genome with columns ['name', 'path', 'taxonomy', 'sp_dir_name']
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

    merged["sp_dir_name"] = (
        merged["taxonomy"].str.split(";").str[-1].str.strip().str.replace(" ", "_")
    )

    matched_df = merged[["name", "path", "taxonomy", "sp_dir_name"]].copy()

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


def apply_species_merging(
    matched_df: pd.DataFrame, species_to_merge_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Reassign taxonomy and sp_dir_name for genomes belonging to merged species.

    Returns a copy of matched_df with updated 'taxonomy', 'sp_dir_name', and an
    'original_taxonomy' column preserving the pre-merge value.
    """
    result = matched_df.copy()
    result["original_taxonomy"] = result["taxonomy"]

    for meta_sp_dir_name, species_to_merge, pangenome_taxonomy in species_to_merge_df[
        ["meta_sp_dir_name", "species_to_merge", "pangenome_taxonomy"]
    ].itertuples(index=False):
        logging.info(
            f"Merging species {species_to_merge} into {meta_sp_dir_name} (metaspecies)"
        )
        for sp in species_to_merge.split(";"):
            filt = result["taxonomy"].str.endswith(f";{sp}")
            result.loc[filt, "taxonomy"] = pangenome_taxonomy
            result.loc[filt, "sp_dir_name"] = meta_sp_dir_name

    logging.info(
        f"After merging, matched_df contains {result['taxonomy'].nunique()} unique taxonomies."
    )
    return result


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
    ].sort_values(
        ["build_pangenome", "input_sp_genome", "species"],
        ascending=[False, False, True],
    )

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
    for (sp_dir_name, taxonomy), group in filtered_df.groupby(
        ["sp_dir_name", "taxonomy"]
    ):
        sp_outdir = outdir / sp_dir_name
        print(">>>>", taxonomy, "||||", sp_dir_name, len(group))
        sp_outdir.mkdir(parents=True, exist_ok=True)
        group[["name", "path"]].sort_values("name").to_csv(
            sp_outdir / "input_genomes.tsv.gz",
            sep="\t",
            index=False,
            header=False,
            compression="gzip",
        )
        with open(sp_outdir / "pangenome_taxonomy.txt", "w") as f:
            f.write(str(taxonomy) + "\n")


def write_metadata_by_species(outdir, genome_metadata_df: pd.DataFrame):
    logging.info(f"Splitting and writing genome metadata by species in {outdir}")
    for sp_dir_name, group in genome_metadata_df.groupby("sp_dir_name"):
        sp_outdir = outdir / sp_dir_name
        sp_outdir.mkdir(parents=True, exist_ok=True)

        group.drop(columns=["sp_dir_name"]).to_csv(
            sp_outdir / "genomes_metadata.tsv.gz",
            sep="\t",
            index=False,
            compression="gzip",
        )


def parse_metadata_file(
    genome_metadata_file: Path, genome_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Parse metadata file and attach sp_dir_name from genome_df.

    genome_df must have columns ['name', 'sp_dir_name'].
    """
    df = pd.read_csv(genome_metadata_file, sep="\t", dtype=str, keep_default_na=False)

    if "genomes" not in df.columns:
        raise KeyError(
            f"genomes column not found in metadata file {genome_metadata_file}"
        )

    df = df[df["genomes"].isin(genome_df["name"])].copy()
    df = df.merge(
        genome_df[["name", "sp_dir_name"]].rename(columns={"name": "genomes"}),
        on="genomes",
        how="left",
    )
    return df


def parse_translation_table_file(
    translation_table_file: Path, genome_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Parse translation table file and associate translation tables with species taxonomy.

    Args:
        translation_table_file: Path to TSV file containing genome accessions and translation tables
        genome_df: DataFrame with columns ['name', 'taxonomy', 'sp_dir_name']

    Returns:
        DataFrame with columns ['name', 'translation_table', 'taxonomy', 'sp_dir_name']
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

    df = df[df["name"].isin(genome_df["name"])].copy()
    df = df.merge(
        genome_df[["name", "taxonomy", "sp_dir_name"]],
        on="name",
        how="left",
    )

    logging.info(
        f"Parsed translation table for {len(df)} genomes across {df['taxonomy'].nunique()} species."
    )
    return df


def compute_pangenome_taxonomy(
    species_to_merge_df: pd.DataFrame, genome_taxonomy_df: pd.DataFrame
) -> pd.DataFrame:
    """For each metaspecies, compute a consensus taxonomy based on the taxonomies of the species to merge."""
    # Map species name (last rank of taxonomy) → full taxonomy string
    tax_unique = genome_taxonomy_df[["taxonomy"]].drop_duplicates().copy()
    tax_unique["species"] = tax_unique["taxonomy"].str.split(";").str[-1].str.strip()
    species_to_tax = tax_unique.set_index("species")["taxonomy"].to_dict()

    def compute_consensus_taxonomy(species_to_merge: str) -> str:
        species_list = [s.strip() for s in species_to_merge.split(";")]
        try:
            full_taxonomies = [species_to_tax[s] for s in species_list]
        except KeyError as e:
            raise ValueError(f"Species {e} not found in taxonomy file") from e

        # Strip GTDB letter suffix from each rank of each taxonomy, then verify consensus
        pangenome_taxonomies = {rm_suffix(tax) for tax in full_taxonomies}
        if len(pangenome_taxonomies) != 1:
            raise ValueError(
                f"Inconsistent pangenome taxonomies for {species_to_merge!r}: {pangenome_taxonomies}"
            )
        return pangenome_taxonomies.pop()

    result = species_to_merge_df.copy()
    result["pangenome_taxonomy"] = result["species_to_merge"].apply(
        compute_consensus_taxonomy
    )
    return result


def write_translation_table_by_species(outdir, translation_df: pd.DataFrame):
    logging.info(f"Splitting and writing translation table by species in {outdir}")
    rows = []
    for (sp_dir_name, taxonomy), group in translation_df.groupby(
        ["sp_dir_name", "taxonomy"]
    ):
        counts = group["translation_table"].value_counts()
        if len(counts) > 1:
            logging.warning(
                f"Species {sp_dir_name} has multiple translation tables: {counts.to_dict()}. The most common one will be used."
            )
        rows.append(f"{sp_dir_name}\t{counts.index[0]}\n")
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


def rm_suffix(species: str) -> str:
    """Remove _[A-Z]+ suffix from species name"""
    return re.sub(r"_[A-Za-z]+$", "", species)


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

    if args.species_to_merge:
        species_to_merge_df = parse_species_to_merge_file(
            args.species_to_merge, args.min_genomes
        )
        species_to_merge_df = compute_pangenome_taxonomy(
            species_to_merge_df, genome_taxonomy_df
        )
        logging.info(
            "Species to merge:\n"
            + species_to_merge_df[
                ["meta_sp_dir_name", "species_to_merge", "pangenome_taxonomy"]
            ].to_csv(sep="\t", index=False)
        )
        matched_df = apply_species_merging(matched_df, species_to_merge_df)

    write_species_summary(
        args.species_summary_file, matched_df, genome_taxonomy_df, args.min_genomes
    )

    # Filter to species suitable for pangenome building
    filtered_df = filter_genomes_for_pangenome(
        matched_df, args.min_genomes, uninformative_species
    )
    write_ppanggolin_input_files(args.outdir, filtered_df)

    if args.genome_metadata:
        metadata_df = parse_metadata_file(args.genome_metadata, filtered_df)
        write_metadata_by_species(args.outdir, metadata_df)

    if args.genome_translation_table:
        translation_df = parse_translation_table_file(
            args.genome_translation_table, filtered_df
        )
        write_translation_table_by_species(args.outdir, translation_df)

if __name__ == "__main__":
    sys.exit(main())
