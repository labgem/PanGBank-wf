#!/usr/bin/env python

import re
import sys
import argparse
import pandas as pd
from pathlib import Path
from collections import defaultdict
from dataclasses import dataclass


def parse_taxonomy(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", names=["genomes", "taxonomy"])


def parse_input_genomes(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", names=["genomes", "path"])


SPLIT_PATTERN = re.compile(r"^s__\S+ \S+_[A-Za-z]+$")

def is_split(species: str) -> bool:
    """Check if a species name is a split"""
    return bool(SPLIT_PATTERN.match(species))

def rm_suffix(species: str) -> str:
    """Remove _[A-Z]+ suffix from species name"""
    return re.sub(r"_[A-Za-z]+$", "", species)

@dataclass
class GTDBMetaSpecies:
    name: str
    splits: dict[str, int]
    genomes: dict[str, list[str]]

    @property
    def size(self) -> int:
        return sum(n for n in self.splits.values())

    def ok(self, min_size: int) -> bool:
        return self.size >= min_size


def find_metaspecies(
    genome_taxonomy: pd.DataFrame, min_size: int
) -> dict[str, GTDBMetaSpecies]:
    """Find all metaspecies. Keeps only metaspecies with more than 'min_size' genomes"""
    species = genome_taxonomy["taxonomy"].str.split(";").str[-1].str.strip()
    metamap = defaultdict(lambda: defaultdict(int))

    for s in species:
        prefix = rm_suffix(s)
        metamap[prefix][s] += 1

    r =  {
        prefix: GTDBMetaSpecies(name=prefix, splits=dict(sc), genomes=defaultdict(list))
        for prefix, sc in metamap.items() if len(sc) > 1
    }

    for _, row in genome_taxonomy.iterrows():
        sp = row["taxonomy"].split(";")[-1].strip()
        meta_sp = rm_suffix(sp)
        if meta_sp in r:
            r[meta_sp].genomes[sp].append(row["genomes"])

    for k in list(r.keys()):
        if sum(len(x) for x in r[k].genomes.values()) < min_size:
            del r[k]

    return r


def filter_pangbank(
    genome_taxonomy: pd.DataFrame,
    metaspecies: dict[str, GTDBMetaSpecies],
    min_size: int,
) -> tuple[pd.DataFrame, list[str]]:
    """Keep all non-split species with more than 'min_size' genomes, i.e. original pangbank filter.
    Keep all metaspecies with more than 'min_size' genomes, i.e. merge candidates.
    """
    species = genome_taxonomy["taxonomy"].str.split(";").str[-1].str.strip()
    species_counts = species.value_counts()

    ok_species = set(species_counts[species_counts >= min_size].index)
    ok_metaspecies = {name for name, ms in metaspecies.items() if ms.ok(min_size)}

    ok_species_no_splits = [s for s in ok_species if rm_suffix(s) not in metaspecies]

    print(
        f"Species with enough genomes (>= {min_size}): {len(ok_species)} total, {len(ok_species_no_splits)} without GTDB splits"
    )
    print(
        f"Metaspecies (split groups with >= {min_size} combined genomes): {len(ok_metaspecies)}/{len(metaspecies)} retained"
    )

    prefixes = species.map(rm_suffix)
    keep = species.isin(ok_species) | prefixes.isin(ok_metaspecies)

    return genome_taxonomy[keep], ok_species_no_splits


def filter_input_genome(genome_df: pd.DataFrame, genomes: set[str]) -> pd.DataFrame:
    filtered = genome_df[genome_df["genomes"].isin(genomes)]
    print(
        f"Restricting taxonomy to input genomes: {len(filtered)}/{len(genome_df)} entries retained "
        f"({len(genome_df) - len(filtered)} not found in taxonomy)"
    )
    return filtered


def filter(
    taxonomy: Path, input_genomes_file: Path, min_genome_count: int
) -> tuple[pd.DataFrame, dict[str, GTDBMetaSpecies], list[str]]:

    genome_taxonomy = parse_taxonomy(taxonomy)
    input_genomes = parse_input_genomes(input_genomes_file)

    print(
        f"Loaded {len(input_genomes)} input genomes and {len(genome_taxonomy)} taxonomy entries"
    )

    genome_taxonomy = filter_input_genome(
        genome_taxonomy, set(input_genomes["genomes"])
    )

    metaspecies = find_metaspecies(genome_taxonomy, min_genome_count)
    genome_taxonomy_filtered, ok_no_splits = filter_pangbank(
        genome_taxonomy, metaspecies, min_genome_count
    )
    print(
        f"Result: {len(ok_no_splits)} non-split pangenomes and {len(metaspecies)} metaspecies to merge"
    )
    return genome_taxonomy_filtered, metaspecies, ok_no_splits


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Find split species using GTDB metadata files.",
        epilog="Example: python find_gtdb_splits.py --taxonomy <path> --input-genomes <path>"
        "--genome-min-completeness <int> "
        "--representative-min-completeness <int> "
        "--min-genome-count <int> --output-directory <path>",
    )

    parser.add_argument(
        "--genome-taxonomy", type=Path, required=True, metavar="FILE", help=""
    )

    parser.add_argument(
        "--input-genomes", type=Path, required=True, metavar="FILE", help=""
    )

    parser.add_argument(
        "--min-genome-count",
        type=int,
        required=True,
        metavar="INT",
        help=""
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
    output = args.output_directory
    output.mkdir(parents=True, exist_ok=True)

    genome_taxonomy_filtered, metaspecies, ok_no_splits = filter(
        args.genome_taxonomy,
        args.input_genomes,
        args.min_genome_count,
    )
    genome_taxonomy_filtered.to_csv(
        output / "taxonomy.filtered.tsv", sep="\t", index=False
    )

    with open(output / "metaspecies.tsv", "w") as fo:
        for k, m in metaspecies.items():
            fo.write(f"{k}\t{';'.join(m.splits)}\n")

    pangenomes = {s: [] for s in ok_no_splits}
    for _, row in genome_taxonomy_filtered.iterrows():
        s = row["taxonomy"].split(";")[-1].strip()
        if s in pangenomes:
            pangenomes[s].append(row["genomes"])

    pangenome_out = output / "pangenomes"
    pangenome_out.mkdir(exist_ok=True)
    for name, genomes in pangenomes.items():
        fn = name.replace(" ", "_") + ".list"
        with open(pangenome_out / fn, "w") as fo:
            fo.writelines(s + "\n" for s in genomes)

    meta_out = output / "meta"
    meta_out.mkdir(exist_ok=True)

    for name, ms in metaspecies.items():
        fn = name.replace(" ", "_") + ".list"
        fns = name.replace(" ", "_") + ".splits"

        with open(meta_out / fn, "w") as fo:
            for split, sps in ms.genomes.items():
                for sp in sps:
                    fo.write(f"{split}\t{sp}\n")

        with open(meta_out / fns, "w") as fo:
            fo.writelines(s + "\n" for s in ms.splits)

if __name__ == "__main__":
    sys.exit(main())
