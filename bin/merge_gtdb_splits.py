#!/usr/bin/env python
import sys
import logging
import argparse
from typing import Any
import pandas as pd
import numpy as np
import networkx as nx
from pathlib import Path
from collections import defaultdict
from itertools import combinations

logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] %(levelname)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)


def build_sketch_to_id(fna_paths_file: Path) -> dict:
    """Map fna path to genome ID from a genome_id<TAB>fna_path file."""
    fna_to_id = {}
    with open(fna_paths_file) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            genome_id, fna_path = line.split("\t", 1)
            fna_to_id[fna_path] = genome_id
    return fna_to_id


def compute_zero_inclusive_median(values: list[float], zero_count: int) -> float:
    """Return the median of positive ANI values plus an implicit number of zeros."""
    total_count = zero_count + len(values)
    if total_count == 0:
        return 0.0

    sorted_values = sorted(values)
    left_index = (total_count - 1) // 2
    right_index = total_count // 2

    def value_at(index: int) -> float:
        if index < zero_count:
            return 0.0
        return sorted_values[index - zero_count]

    return (value_at(left_index) + value_at(right_index)) / 2


def build_species_ani_df(
    skani_dist_path: Path,
    accession_to_species: dict[str, str],
    species_to_accessions: dict[str, list[str]],
    fna_paths_file: Path,
    af_threshold: float,
    matrix_ani_stat: str = "mean",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build a species-level mean ANI matrix by streaming a skani dist TSV file.

    All possible species pairs are inferred from ``species_to_accessions``.
    Pairs absent from the file are assumed to be too divergent and treated as
    ANI=0, AF=0. Pairs present but below ``af_threshold`` are also treated as
    ANI=0, AF=0. In both cases the pair still contributes to the denominator
    when computing per-species-pair means, so ``mean_ani`` and ``mean_af`` are
    always averaged over ``theoretical_pair_count``.

    Args:
        skani_dist_path: Path to the tab-separated skani dist output file.
            May be empty (header only), in which case all pairs are treated
            as ANI=0.
        accession_to_species: Mapping from genome accession to species label.
        species_to_accessions: Mapping from species label to its list of
            accessions. Used to enumerate all theoretical species pairs.
        fna_paths_file: File listing sketch paths used to build the accession
            lookup.
        af_threshold: Minimum alignment fraction for a pair to contribute a
            non-zero ANI/AF value. Pairs below this threshold are treated as
            ANI=0, AF=0 (same as missing pairs).
        matrix_ani_stat: Summary statistic stored in species_ani_df. Must be
            either "mean" or "median".

    Returns:
        species_ani_df: Square DataFrame (species × species) of ANI summary
            values averaged over all theoretical pairs (missing and low-AF
            pairs contribute 0 for both mean and median calculations).
            Diagonal and cells with no observed pairs are 0.
        pair_summary_df: Long-form DataFrame with one row per species pair,
            containing mean ANI, median ANI, mean AF, and per-category pair
            counts for QC.
            Columns:
                species_A, species_B        : species labels (sorted)
                mean_ani                    : mean ANI over theoretical_pair_count
                median_ani                  : median ANI over theoretical_pair_count
                mean_af                     : mean AF over theoretical_pair_count
                theoretical_pair_count      : total possible genome pairs
                observed_count              : pairs present in the skani file
                passing_af_count            : observed pairs that met af_threshold
                below_af_count              : observed pairs below af_threshold
                missing_count               : theoretical pairs absent from file
    """

    if matrix_ani_stat not in {"mean", "median"}:
        raise ValueError(
            f"Unsupported ANI matrix summary statistic: {matrix_ani_stat!r}. "
            "Expected 'mean' or 'median'."
        )

    log.info(f"Building sketch-path-to-accession map from: {fna_paths_file}")
    sketch_path_to_accession = build_sketch_to_id(fna_paths_file)
    log.info(
        f"Mapped {len(sketch_path_to_accession)} sketch paths to genome accessions"
    )

    # Pre-populate all theoretical species pairs from species_to_accessions so
    # that pairs entirely absent from the skani file are still represented.
    all_species = sorted(species_to_accessions.keys())
    species_pair_stats: dict[tuple[str, str], dict[str, Any]] = {}
    for sp_a, sp_b in combinations(all_species, 2):
        species_pair: tuple[str, str] = tuple(sorted((sp_a, sp_b)))  # type: ignore[assignment]
        species_pair_stats[species_pair] = {
            "ani_sum": 0.0,
            "af_sum": 0.0,
            "passing_ani_values": [],
            "observed_count": 0,
            "passing_af_count": 0,
            "below_af_count": 0,
        }

    log.info(
        f"Initialized {len(species_pair_stats)} theoretical species pairs "
        f"across {len(all_species)} species"
    )

    log.info(f"Streaming skani dist file: {skani_dist_path}")
    n_lines = 0
    n_skipped_nan = 0
    n_below_af = 0

    with open(skani_dist_path) as f:
        header = f.readline().rstrip("\n").split("\t")
        ref_col = header.index("Ref_file")
        query_col = header.index("Query_file")
        ani_col = header.index("ANI")
        query_af_col = header.index("Align_fraction_query")
        ref_af_col = header.index("Align_fraction_ref")

        for line in f:
            n_lines += 1
            parts = line.rstrip("\n").split("\t")

            ref_accession = sketch_path_to_accession.get(parts[ref_col])
            query_accession = sketch_path_to_accession.get(parts[query_col])
            if ref_accession is None or query_accession is None:
                raise ValueError(
                    f"Could not map sketch path to accession in line: {line!r}"
                )

            ani = float(parts[ani_col])
            if np.isnan(ani):
                n_skipped_nan += 1
                continue

            query_af = float(parts[query_af_col])
            ref_af = float(parts[ref_af_col])
            pair_af = max(query_af, ref_af)

            sp_ref = accession_to_species.get(ref_accession)
            sp_query = accession_to_species.get(query_accession)
            if sp_ref is None or sp_query is None:
                raise ValueError(
                    f"Could not map accession '{ref_accession}' or '{query_accession}' "
                    f"to a species in line: {line!r}"
                )

            species_pair = tuple(sorted((sp_ref, sp_query)))  # type: ignore[assignment]

            if species_pair not in species_pair_stats:
                raise ValueError(
                    f"Species pair {species_pair} found in skani output but not in theoretical pairs. "
                    f"This likely means accession '{ref_accession}' or '{query_accession}' "
                    f"is missing from species_to_accessions. Line: {line!r}"
                )

            species_pair_stats[species_pair]["observed_count"] += 1

            if pair_af < af_threshold:
                n_below_af += 1
                species_pair_stats[species_pair]["below_af_count"] += 1
            else:
                species_pair_stats[species_pair]["ani_sum"] += ani
                species_pair_stats[species_pair]["af_sum"] += pair_af
                species_pair_stats[species_pair]["passing_ani_values"].append(ani)
                species_pair_stats[species_pair]["passing_af_count"] += 1

    log.info(
        f"Finished streaming {n_lines} data lines — "
        f"{n_skipped_nan} skipped (NaN ANI), "
        f"{n_below_af} treated as zero-ANI (AF < {af_threshold})"
    )

    sorted_species = all_species
    log.info(
        f"Building {len(sorted_species)} × {len(sorted_species)} species ANI matrix"
    )

    species_ani_df = pd.DataFrame(
        index=sorted_species, columns=sorted_species, dtype=float
    )

    pair_summary_records: list[dict[str, Any]] = []

    for species_pair, stats in species_pair_stats.items():
        sp_a, sp_b = species_pair
        theoretical_pair_count = len(species_to_accessions[sp_a]) * len(
            species_to_accessions[sp_b]
        )

        # Missing pairs (not in file) and below-AF pairs both contribute 0 to
        # sums, so dividing by theoretical_pair_count captures all three cases:
        # passing pairs (full ANI/AF), below-AF pairs (0), missing pairs (0).
        mean_ani = stats["ani_sum"] / theoretical_pair_count
        mean_af = stats["af_sum"] / theoretical_pair_count
        zero_ani_count = theoretical_pair_count - stats["passing_af_count"]
        median_ani = compute_zero_inclusive_median(
            stats["passing_ani_values"], zero_ani_count
        )
        missing_count = theoretical_pair_count - stats["observed_count"]

        pair_summary_records.append(
            {
                "species_A": sp_a,
                "species_B": sp_b,
                "mean_ani": mean_ani,
                "median_ani": median_ani,
                "mean_af": mean_af,
                "theoretical_pair_count": theoretical_pair_count,
                "observed_count": stats["observed_count"],
                "passing_af_count": stats["passing_af_count"],
                "below_af_count": stats["below_af_count"],
                "missing_count": missing_count,
            }
        )

        matrix_ani_value = mean_ani if matrix_ani_stat == "mean" else median_ani
        species_ani_df.loc[sp_a, sp_b] = matrix_ani_value
        species_ani_df.loc[sp_b, sp_a] = matrix_ani_value

    species_ani_df = species_ani_df.fillna(0)
    pair_summary_df = pd.DataFrame(pair_summary_records)

    total_missing = (
        pair_summary_df["missing_count"].sum() if not pair_summary_df.empty else 0
    )
    log.info(
        f"Species ANI matrix built — "
        f"{total_missing:.0f} missing genome pairs across all species pairs "
        f"(treated as ANI=0)"
    )

    return species_ani_df, pair_summary_df


def construct_graph(df: pd.DataFrame, ani_threshold: float):
    log.info(f"Building ANI graph with threshold {ani_threshold:.4f}")
    G = nx.Graph()
    for sp in df.index:
        G.add_node(sp)

    for i, sp_a in enumerate(df.index):
        for j, sp_b in enumerate(df.index):
            if j <= i:
                continue
            v = df.loc[sp_a, sp_b]
            if pd.notna(v) and v >= ani_threshold:
                G.add_edge(sp_a, sp_b)

    log.info(f"Graph has {G.number_of_nodes()} nodes and {G.number_of_edges()} edges")
    return G


def find_clusters(G: nx.Graph):
    clusters = list(nx.connected_components(G))
    log.info(f"Found {len(clusters)} connected components (clusters)")
    return clusters


def parse_gtdb_metadata(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")


def parse_genome_list(path: Path) -> tuple[dict[str, str], dict[str, list[str]]]:
    """Parse a two-column TSV (species<TAB>accession) and return
    accession→species and species→[accessions] mappings."""
    log.info(f"Parsing genome list: {path}")
    acc_to_species = {}
    species_to_acc = defaultdict(list)
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            species, acc = line.split("\t")
            acc_to_species[acc] = species
            species_to_acc[species].append(acc)
    log.info(
        f"Loaded {len(acc_to_species)} genomes across {len(species_to_acc)} species"
    )
    return acc_to_species, species_to_acc


def write_clusters(
    clusters, species_to_acc: dict[str, list[str]], prefix: str
) -> tuple[int, int, int]:
    """Write <prefix>.clusters and <prefix>.genomes.clusters files.
    Returns (n_merged_splits, n_singleton_splits, n_merged_clusters)."""
    n_cluster_index = 0
    n_merged_splits = 0
    n_singleton_splits = 0
    n_merged_clusters = 0
    with open(f"{prefix}.clusters", "w") as f_clusters, open(
        f"{prefix}.genomes.clusters", "w"
    ) as f_genomes:
        for cluster in clusters:
            cluster_genomes = [
                acc for split in cluster for acc in species_to_acc[split]
            ]

            if len(cluster) == 1:
                cluster_name = list(cluster)[0].replace(" ", "_")
                n_singleton_splits += 1
            else:
                cluster_name = f"{prefix}_{n_cluster_index}"
                n_cluster_index += 1
                n_merged_splits += len(cluster)
                n_merged_clusters += 1

            f_clusters.write(
                f"{cluster_name}\t{';'.join(cluster)}\t{len(cluster_genomes)}\n"
            )
            f_genomes.write(
                f"{cluster_name}\t{';'.join(cluster_genomes)}\t{len(cluster_genomes)}\n"
            )

    return n_merged_splits, n_singleton_splits, n_merged_clusters


def write_merge_summary(
    prefix: str, n_merged_splits: int, n_singleton_splits: int, n_merged_clusters: int
):
    """Write a one-row TSV summary of the merging statistics."""
    with open(f"{prefix}.merge_summary.tsv", "w") as f:
        f.write("metaspecies\tmerged_splits\tsingleton_splits\tmerged_clusters\n")
        f.write(
            f"{prefix}\t{n_merged_splits}\t{n_singleton_splits}\t{n_merged_clusters}\n"
        )


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Merge GTDB split species",
        epilog="Example: python merge_gtdb_splits.py --genome-list <path> --skani-triangle <path> --ani_threshold <float> --af_threshold <float> --prefix <str>",
    )

    parser.add_argument(
        "--genome-list", type=Path, required=True, metavar="FILE", help=""
    )

    parser.add_argument(
        "--skani-dist", type=Path, required=True, metavar="FILE", help=""
    )

    parser.add_argument(
        "--genome-fna-paths",
        type=Path,
        required=True,
        metavar="FILE",
        help="Two-column TSV: genome_id<TAB>fna_path; used to map sketch paths back to genome IDs",
    )

    parser.add_argument(
        "--ani_threshold", type=float, required=True, metavar="FLOAT", help=""
    )

    parser.add_argument(
        "--af_threshold", type=float, required=True, metavar="FLOAT", help=""
    )

    parser.add_argument(
        "--species-ani-stat",
        type=str,
        default="mean",
        choices=("mean", "median"),
        metavar="{mean,median}",
        help="Summary statistic used to populate the species ANI matrix.",
    )

    parser.add_argument("--prefix", type=str, required=True, metavar="STR", help="")

    return parser.parse_args(argv)


def main():
    args = parse_args()

    acc_to_species, species_to_acc = parse_genome_list(args.genome_list)

    species_ani_df, pair_info_df = build_species_ani_df(
        args.skani_dist,
        acc_to_species,
        species_to_acc,
        args.genome_fna_paths,
        args.af_threshold,
        args.species_ani_stat,
    )

    log.info(
        f"Saving species ANI matrix ({args.species_ani_stat}) to "
        f"{args.prefix}.species_ani.tsv"
    )
    species_ani_df.to_csv(f"{args.prefix}.species_ani.tsv", sep="\t")

    log.info(f"Saving species pair summary to {args.prefix}.species_pair_summary.tsv")
    columns = pair_info_df.columns
    pair_info_df["main_species"] = args.prefix
    pair_info_df[["main_species"] + list(columns)].to_csv(
        f"{args.prefix}.species_pair_summary.tsv", sep="\t", index=False
    )

    species_ani_graph = construct_graph(species_ani_df, float(args.ani_threshold))
    clusters = find_clusters(species_ani_graph)

    n_merged_splits, n_singleton_splits, n_merged_clusters = write_clusters(
        clusters, species_to_acc, args.prefix
    )
    log.info(
        f"Results: {n_merged_clusters} merged clusters ({n_merged_splits} splits merged), {n_singleton_splits} singleton splits"
    )
    write_merge_summary(
        args.prefix, n_merged_splits, n_singleton_splits, n_merged_clusters
    )
    log.info(f"Done. Output written with prefix '{args.prefix}'")


if __name__ == "__main__":
    sys.exit(main())
