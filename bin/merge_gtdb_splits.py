#!/usr/bin/env python
import sys
import logging
import argparse
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
    fna_paths_file: Path,
) -> pd.DataFrame:
    """Stream the skani dist TSV and directly accumulate per-species-pair ANI
    statistics (sum + count). No genome-level matrix or ANI list is kept in
    memory.  Returns a species x species DataFrame of mean ANI values."""
    log.info(f"Building sketch-to-genome-ID map from: {fna_paths_file}")
    fna_to_id = build_sketch_to_id(fna_paths_file)
    log.info(f"Mapped {len(fna_to_id)} sketch paths to genome IDs")

    # ani_sum[pair] and ani_count[pair] accumulate the running mean components
    ani_sum: dict[tuple, float] = defaultdict(float)
    ani_count: dict[tuple, int] = defaultdict(int)
    species_set: set[str] = set()

    log.info(f"Streaming skani dist file: {skani_dist_path}")
    with open(skani_dist_path) as f:
        header = f.readline().rstrip("\n").split("\t")
        ref_col = header.index("Ref_file")
        query_col = header.index("Query_file")
        ani_col = header.index("ANI")

        for line in f:
            parts = line.rstrip("\n").split("\t")
            ref = fna_to_id.get(parts[ref_col])
            query = fna_to_id.get(parts[query_col])
            if ref is None or query is None:
                raise ValueError(
                    f"Could not map sketch path to genome ID in line: {line}"
                )

            ani = float(parts[ani_col])

            if np.isnan(ani):
                continue

            sp_ref = accession_to_species.get(ref)
            sp_query = accession_to_species.get(query)

            if sp_ref is None or sp_query is None:
                raise ValueError(
                    f"Could not map genome ID {ref} or {query} to species in line: {line}"
                )

            species_pair = tuple(sorted((sp_ref, sp_query)))
            ani_sum[species_pair] += ani
            ani_count[species_pair] += 1
            species_set.update([sp_ref, sp_query])

    n_pairs = len(ani_sum)
    log.info(
        f"Processed {n_pairs} inter-species pairs across {len(species_set)} species"
    )
    log.info(
        f"Building species ANI DataFrame ({len(species_set)} x {len(species_set)})"
    )
    sorted_species = sorted(species_set)
    species_ani_df = pd.DataFrame(
        index=sorted_species, columns=sorted_species, dtype=float
    )

    for species_pair, total in ani_sum.items():
        sp_a, sp_b = species_pair
        mean_ani = total / ani_count[species_pair]
        species_ani_df.loc[sp_a, sp_b] = mean_ani
        species_ani_df.loc[sp_b, sp_a] = mean_ani

    species_ani_df = species_ani_df.fillna(0)

    return species_ani_df


def construct_graph(df: pd.DataFrame, threshold: float):
    log.info(f"Building ANI graph with threshold {threshold:.4f}")
    G = nx.Graph()
    for sp in df.index:
        G.add_node(sp)

    for i, sp_a in enumerate(df.index):
        for j, sp_b in enumerate(df.index):
            if j <= i:
                continue
            v = df.loc[sp_a, sp_b]
            if pd.notna(v) and v >= threshold:
                G.add_edge(sp_a, sp_b)

    log.info(f"Graph has {G.number_of_nodes()} nodes and {G.number_of_edges()} edges")
    return G

def find_clusters(G: nx.Graph):
    clusters = list(nx.connected_components(G))
    log.info(f"Found {len(clusters)} connected components (clusters)")
    return clusters

def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Merge GTDB split species",
        epilog="Example: python merge_gtdb_splits.py --genome-list <path> --skani-triangle <path> --threshold <float> --prefix <str>"
    )

    parser.add_argument(
        "--genome-list",
        type=Path,
        required=True,
        metavar="FILE",
        help=""
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
        "--threshold",
        type=float,
        required=True,
        metavar="FLOAT",
        help=""
    )

    parser.add_argument(
        "--prefix",
        type=str,
        required=True,
        metavar="STR",
        help=""
    )
    return parser.parse_args(argv)

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


def main():
    args = parse_args()

    acc_to_species, species_to_acc = parse_genome_list(args.genome_list)

    species_ani_df = build_species_ani_df(
        args.skani_dist, acc_to_species, args.genome_fna_paths
    )

    log.info(f"Saving species ANI matrix to {args.prefix}.species_ani.tsv")
    species_ani_df.to_csv(f"{args.prefix}.species_ani.tsv", sep="\t")

    species_ani_graph = construct_graph(species_ani_df, float(args.threshold))
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
