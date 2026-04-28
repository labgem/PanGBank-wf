import re
import sys
import argparse
import pandas as pd
import numpy as np
import networkx as nx
from pathlib import Path
from collections import defaultdict
from dataclasses import dataclass

def parse_skani_triangle(path):
    with open(path) as f:
        n = int(f.readline().strip())
        labels = []
        vals = np.full((n, n), np.nan)

        for i in range(n):
            parts = f.readline().strip().split("\t")

            label = parts[0]
            values = parts[1:]

            labels.append(label)

            for j, v in enumerate(values):
                vals[i, j] = float(v)

    return labels, vals

def construct_df(labels: list[str], matrix, accession_to_species: dict[str, str]):
    n = len(labels)
    species = set()
    pairs = defaultdict(list)

    for i in range(n):
        for j in range(i + 1, n):
            a, b = labels[i], labels[j]
            sp_a, sp_b = accession_to_species[a], accession_to_species[b]
            species.update([sp_a, sp_b])
            v, k = matrix[j, i], tuple(sorted((sp_a, sp_b)))
            pairs[k].append(v)

    s_species = sorted(species)

    df = pd.DataFrame(index=s_species, columns=s_species, dtype=float)

    for (sp_a, sp_b), vals in pairs.items():
        if len(vals):
            m = np.mean(vals)
            df.loc[sp_a, sp_b] = m
            df.loc[sp_b, sp_a] = m
        else:
            df.loc[sp_a, sp_b] = np.nan
            df.loc[sp_b, sp_a] = np.nan

    for sp in s_species:
        intra = pairs.get((sp, sp), [])
        df.loc[sp, sp] = np.mean(intra) if intra else np.nan

    df = df.fillna(0)

    return df

def construct_graph(df: pd.DataFrame, threshold: float):
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

    return G

def find_clusters(G: nx.Graph):
    return list(nx.connected_components(G))

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
        "--skani-triangle",
        type=Path,
        required=True,
        metavar="FILE",
        help=""
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

def main():
    args = parse_args()

    labels, matrix = parse_skani_triangle(args.skani_triangle)

    acc_to_species = {}
    species_to_acc = defaultdict(list)

    with open(args.genome_list) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            species, acc = line.strip().split("\t")
            acc_to_species[acc] = species
            species_to_acc[species].append(acc)

    df = construct_df(labels, matrix, acc_to_species)
    G = construct_graph(df, float(args.threshold))
    clusters = find_clusters(G)

    with open(f"{args.prefix}.clusters", "w") as f, open(f"{args.prefix}.genomes.clusters", "w") as f2:
        for i, cluster in enumerate(clusters):
            f.write(f"{args.prefix}_{i}\t")
            s = ",".join(cluster)
            f.write(s + "\n")

            f2.write(f"{args.prefix}_{i}\t")
            all_genomes = []
            for split in cluster:
                all_genomes.extend(species_to_acc[split])
            f2.write(",".join(all_genomes))

if __name__ == "__main__":
    sys.exit(main())



