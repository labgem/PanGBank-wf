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


def build_species_ani_df(
    skani_dist_path: Path,
    accession_to_species: dict[str, str],
    species_to_accessions: dict[str, list[str]],
    fna_paths_file: Path,
    af_threshold: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build a species-level mean ANI matrix by streaming a skani dist TSV file.

    Pairs absent from the file are assumed to be too divergent and treated as
    ANI=0, AF=0. Pairs present but below ``af_threshold`` are also treated as
    ANI=0, AF=0. In both cases the pair still contributes to the denominator
    when computing per-species-pair means, so ``mean_ani`` and ``mean_af`` are
    always averaged over ``theoretical_pair_count``.

    Args:
        skani_dist_path: Path to the tab-separated skani dist output file.
        accession_to_species: Mapping from genome accession to species label.
        species_to_accessions: Mapping from species label to its list of accessions.
        fna_paths_file: File listing sketch paths used to build the accession lookup.
        af_threshold: Minimum alignment fraction for a pair to contribute a
            non-zero ANI/AF value. Pairs below this threshold are treated as
            ANI=0, AF=0 (same as missing pairs).

    Returns:
        species_ani_df: Square DataFrame (species × species) of mean ANI values
            averaged over all theoretical pairs (missing and low-AF pairs
            contribute 0). Diagonal and cells with no observed pairs are 0.
        pair_summary_df: Long-form DataFrame with one row per species pair,
            containing mean ANI, mean AF, and per-category pair counts for QC.
            Columns:
                species_A, species_B        : species labels (sorted)
                mean_ani                    : mean ANI over theoretical_pair_count
                mean_af                     : mean AF over theoretical_pair_count
                theoretical_pair_count      : total possible genome pairs
                observed_count              : pairs present in the skani file
                passing_af_count            : observed pairs that met af_threshold
                below_af_count              : observed pairs below af_threshold
                missing_count               : theoretical pairs absent from file
    """

    log.info(f"Building sketch-path-to-accession map from: {fna_paths_file}")
    sketch_path_to_accession = build_sketch_to_id(fna_paths_file)
    log.info(
        f"Mapped {len(sketch_path_to_accession)} sketch paths to genome accessions"
    )

    # Keyed by sorted (species_A, species_B) tuple to ensure symmetry.
    species_pair_stats: dict[tuple[str, str], dict[str, Any]] = {}
    observed_species: set[str] = set()

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

            species_pair: tuple[str, str] = tuple(sorted((sp_ref, sp_query)))  # type: ignore[assignment]

            if species_pair not in species_pair_stats:
                species_pair_stats[species_pair] = {
                    "ani_sum": 0.0,
                    "af_sum": 0.0,
                    "observed_count": 0,
                    "passing_af_count": 0,
                    "below_af_count": 0,
                }

            species_pair_stats[species_pair]["observed_count"] += 1
            observed_species.update([sp_ref, sp_query])

            if pair_af < af_threshold:
                # Treated as ANI=0, AF=0 — contributes nothing to sums but
                # is counted in the denominator via theoretical_pair_count.
                n_below_af += 1
                species_pair_stats[species_pair]["below_af_count"] += 1
            else:
                species_pair_stats[species_pair]["ani_sum"] += ani
                species_pair_stats[species_pair]["af_sum"] += pair_af
                species_pair_stats[species_pair]["passing_af_count"] += 1

    log.info(
        f"Finished streaming {n_lines} data lines — "
        f"{n_skipped_nan} skipped (NaN ANI), "
        f"{n_below_af} treated as zero-ANI (AF < {af_threshold}), "
        f"{len(species_pair_stats)} unique species pairs observed"
    )

    sorted_species = sorted(observed_species)
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
        missing_count = theoretical_pair_count - stats["observed_count"]

        pair_summary_records.append(
            {
                "species_A": sp_a,
                "species_B": sp_b,
                "mean_ani": mean_ani,
                "mean_af": mean_af,
                "theoretical_pair_count": theoretical_pair_count,
                "observed_count": stats["observed_count"],
                "passing_af_count": stats["passing_af_count"],
                "below_af_count": stats["below_af_count"],
                "missing_count": missing_count,
            }
        )

        species_ani_df.loc[sp_a, sp_b] = mean_ani
        species_ani_df.loc[sp_b, sp_a] = mean_ani

    species_ani_df = species_ani_df.fillna(0)
    pair_summary_df = pd.DataFrame(pair_summary_records)

    log.info(
        f"Species ANI matrix built — "
        f"{pair_summary_df['missing_count'].sum():.0f} missing pairs across all species pairs "
        f"(treated as ANI=0)"
    )

    return species_ani_df, pair_summary_df


def build_species_ani_df_old(
    skani_dist_path: Path,
    accession_to_species: dict[str, str],
    species_to_accessions: dict[str, list[str]],
    fna_paths_file: Path,
    af_threshold: float,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Stream the skani dist TSV and directly accumulate per-species-pair ANI
    statistics (sum + count). No genome-level matrix or ANI list is kept in
    memory.  Returns a species x species DataFrame of mean ANI values."""
    log.info(f"Building sketch-to-genome-ID map from: {fna_paths_file}")
    fna_to_id = build_sketch_to_id(fna_paths_file)
    log.info(f"Mapped {len(fna_to_id)} sketch paths to genome IDs")

    species_pair_info: dict[tuple[str, ...], dict[str, Any]] = {}

    species_set: set[str] = set()

    log.info(f"Streaming skani dist file: {skani_dist_path}")

    with open(skani_dist_path) as f:
        header = f.readline().rstrip("\n").split("\t")
        ref_col = header.index("Ref_file")
        query_col = header.index("Query_file")
        ani_col = header.index("ANI")
        query_af_col = header.index("Align_fraction_query")
        ref_af_col = header.index("Align_fraction_ref")

        for line in f:
            parts = line.rstrip("\n").split("\t")
            ref = fna_to_id.get(parts[ref_col])
            query = fna_to_id.get(parts[query_col])
            if ref is None or query is None:
                raise ValueError(
                    f"Could not map sketch path to genome ID in line: {line}"
                )

            ani = float(parts[ani_col])
            query_af = float(parts[query_af_col])
            ref_af = float(parts[ref_af_col])

            af = max(query_af, ref_af)

            if np.isnan(ani):
                continue

            sp_ref = accession_to_species.get(ref)
            sp_query = accession_to_species.get(query)

            if sp_ref is None or sp_query is None:
                raise ValueError(
                    f"Could not map genome ID {ref} or {query} to species in line: {line}"
                )

            species_pair = tuple(sorted((sp_ref, sp_query)))

            if species_pair not in species_pair_info:

                species_pair_info[species_pair] = {
                    "species_pair": species_pair,
                    "ani_sum": 0.0,
                    "ani_count": 0,
                    "low_af_count_filtered": 0,
                    "af_sum": 0.0,
                    "af_count": 0,
                }

            if af < af_threshold:
                species_pair_info[species_pair]["low_af_count_filtered"] += 1
            else:
                species_pair_info[species_pair]["ani_sum"] += ani
                species_pair_info[species_pair]["ani_count"] += 1
                species_pair_info[species_pair]["af_sum"] += af
                species_pair_info[species_pair]["af_count"] += 1

            # ani_sum[species_pair] += ani
            # ani_count[species_pair] += 1
            species_set.update([sp_ref, sp_query])

    n_pairs = len(species_pair_info)
    log.info(
        f"Processed {n_pairs} inter-species pairs across {len(species_set)} species"
    )
    # log.info(
    #     f"Building species ANI DataFrame ({len(species_set)} x {len(species_set)})"
    # )

    sorted_species = sorted(species_set)
    species_ani_df = pd.DataFrame(
        index=sorted_species, columns=sorted_species, dtype=float
    )
    pair_info_list = []
    for species_pair, info in species_pair_info.items():
        sp_a, sp_b = species_pair
        theoretical_pair_count = len(species_to_accessions[sp_a]) * len(
            species_to_accessions[sp_b]
        )

        mean_ani = info["ani_sum"] / theoretical_pair_count
        mean_af = info["af_sum"] / theoretical_pair_count
        summary_info = {
            "species_A": sp_a,
            "species_B": sp_b,
            "mean_ani": mean_ani,
            "mean_af": mean_af,
            "theoretical_pair_count": theoretical_pair_count,
            "low_af_count_filtered": info.get("low_af_count_filtered", 0),
            "ani_count": info.get("ani_count", 0),
        }

        pair_info_list.append(summary_info)

        species_ani_df.loc[sp_a, sp_b] = mean_ani
        species_ani_df.loc[sp_b, sp_a] = mean_ani

    species_ani_df = species_ani_df.fillna(0)

    pair_info_df = pd.DataFrame(pair_info_list)

    return species_ani_df, pair_info_df


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
    )

    log.info(f"Saving species ANI matrix to {args.prefix}.species_ani.tsv")
    species_ani_df.to_csv(f"{args.prefix}.species_ani.tsv", sep="\t")

    log.info(f"Saving species pair summary to {args.prefix}.species_pair_summary.tsv")
    pair_info_df.to_csv(
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
