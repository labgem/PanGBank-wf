#!/usr/bin/env python
import argparse
import pandas as pd
from pathlib import Path

def parse_gtdb_metadata(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")

def parse_input_genomes(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", names=["genomes", "path"])


def filter_representative(metadata: pd.DataFrame, min_checkm: float) -> pd.DataFrame:
    """Apply checkm/checkm2 completeness filtering on representative genomes"""
    rm = set()
    nr = 0
    for _, row in metadata.iterrows():
        acc = row["genomes"]
        if acc == row["gtdb_genome_representative"]:
            nr += 1
            cm, cm2 = row["checkm_completeness"], row["checkm2_completeness"]
            if max(cm, cm2) < min_checkm:
                rm.add(acc)

    filtered = metadata[
        ~metadata["genomes"].isin(rm) &
        ~metadata["gtdb_genome_representative"].isin(rm)
    ]

    print(f"Representative filtering (min completeness: {min_checkm}%): {len(rm)}/{nr} bad representatives removed, discarding {len(metadata)-len(filtered)}/{len(metadata)} genomes")

    return filtered

def filter_genome(metadata: pd.DataFrame, min_checkm: float) -> pd.DataFrame:
    """Apply checkm/checkm2 completeness filtering on all genomes"""
    filtered = metadata[
        metadata[["checkm_completeness", "checkm2_completeness"]].max(axis=1) >= min_checkm
    ]
    print(f"Genome filtering (min completeness: {min_checkm}%): {len(metadata) - len(filtered)}/{len(metadata)} genomes removed")
    return filtered

def filter_input_genome(metadata: pd.DataFrame, genomes: set[str], label: str = "Genome input filtering") -> pd.DataFrame:
    filtered = metadata[metadata["genomes"].isin(genomes)]
    print(f"{label}: {len(filtered)}/{len(metadata)} genomes retained")
    return filtered

def check_input_genomes(input_genomes: pd.DataFrame) -> None:
    """Check that genome names and paths are unique in the input file."""
    dup_names = input_genomes["genomes"][input_genomes["genomes"].duplicated()].unique()
    if len(dup_names) > 0:
        raise ValueError(
            f"Duplicate genome names found in input file: {', '.join(dup_names)}"
        )

    dup_paths = input_genomes["path"][input_genomes["path"].duplicated()].unique()
    if len(dup_paths) > 0:
        raise ValueError(
            f"Duplicate genome paths found in input file: {', '.join(dup_paths)}"
        )


def filter_genomes(
    metadata_path: Path,
    input_genomes_file: Path,
    min_checkm_completeness_repr: float,
    min_checkm_completeness: float,
) -> tuple[pd.DataFrame, pd.DataFrame, dict]:

    input_genomes = parse_input_genomes(input_genomes_file)
    check_input_genomes(input_genomes)
    metadata = parse_gtdb_metadata(metadata_path)

    n_input = len(input_genomes)
    print(f"Starting quality filtering: {n_input} input genomes, {len(metadata)} entries in GTDB metadata")

    metadata = filter_input_genome(metadata, set(input_genomes["genomes"]), label="Metadata restricted to input genomes")
    n_in_metadata = len(metadata)
    n_not_in_metadata = n_input - n_in_metadata

    if n_not_in_metadata > 0:
        missing = set(input_genomes["genomes"]) - set(metadata["genomes"])
        raise ValueError(
            f"{n_not_in_metadata} input genome(s) not found in metadata: {', '.join(sorted(missing))}"
        )

    if "gtdb_genome_representative" in metadata.columns:
        metadata = filter_representative(metadata, min_checkm_completeness_repr)
    else:
        print(
            "Warning: no 'gtdb_genome_representative' column found in metadata, skipping representative filtering"
        )
    n_after_repr = len(metadata)
    n_representative_filtered = n_in_metadata - n_after_repr

    metadata = filter_genome(metadata, min_checkm_completeness)
    n_after_genome = len(metadata)
    n_genome_quality_filtered = n_after_repr - n_after_genome

    input_genomes = filter_input_genome(input_genomes, set(metadata["genomes"]), label="Input genomes matching quality-filtered metadata")
    n_retained = len(input_genomes)

    print(f"Quality filtering complete: {n_retained}/{n_input} genomes retained")

    stats = {
        "input_genomes": n_input,
        "retained_genomes": n_retained,
        "representative_filtered": n_representative_filtered,
        "genome_quality_filtered": n_genome_quality_filtered,
    }

    return metadata, input_genomes, stats


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Filter genomes based on quality using GTDB metadata files.",
        epilog="Example: python filter_genome_on_quality.py --metadata-file <path> --used-genomes <path>"
        "--genome-min-completeness <int> "
        "--representative-min-completeness <int> "
        "--min-genome-count <int> --output-directory <path>",
    )

    parser.add_argument(
        "--metadata-file",
        type=Path,
        required=True,
        metavar="FILE",
        help=""
    )

    parser.add_argument(
        "--input-genomes",
        type=Path,
        required=True,
        metavar="FILE",
        help=""
    )

    parser.add_argument(
        "--genome-min-completeness", type=int, required=True, metavar="INT", help=""
    )

    parser.add_argument(
        "--representative-min-completeness",
        type=int,
        required=True,
        metavar="INT",
        help="",
    )

    parser.add_argument(
        "--output-directory",
        type=Path,
        required=True,
        metavar="DIRECTORY",
        help=""
    )

    return parser.parse_args(argv)

def main():
    args = parse_args()
    outdir = args.output_directory
    outdir.mkdir(parents=True, exist_ok=True)

    filtered_metadata, filtered_input_genomes, stats = filter_genomes(
        args.metadata_file,
        args.input_genomes,
        args.representative_min_completeness,
        args.genome_min_completeness,
    )
    filtered_metadata.to_csv(outdir / "genome_metadata.filtered.tsv", sep="\t", index=False)
    filtered_input_genomes.to_csv(outdir / "input_genomes.filtered.tsv", sep="\t", index=False, header=False)

    with open(outdir / "genome_quality_filtering_summary.tsv", "w") as f:
        f.write(
            "sample\tretained_genomes\trepresentative_filtered\tgenome_quality_filtered\n"
        )
        f.write(
            f"genome_quality_filtering\t{stats['retained_genomes']}\t"
            f"{stats['representative_filtered']}\t{stats['genome_quality_filtered']}\n"
        )


if __name__ == "__main__":
    main()
