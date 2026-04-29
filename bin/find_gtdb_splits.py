#!/usr/bin/env python

import re
import sys
import argparse
import pandas as pd
from pathlib import Path
from collections import defaultdict
from dataclasses import dataclass

def parse_gtdb_metadata(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")

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

    print(f"Representative filtering: {len(rm)}/{nr} bad representatives, removing {len(metadata)-len(filtered)}/{len(metadata)} genomes")

    return filtered

def filter_genome(metadata: pd.DataFrame, min_checkm: float) -> pd.DataFrame:
    """Apply checkm/checkm2 completeness filtering on all genomes"""
    filtered = metadata[
        metadata[["checkm_completeness", "checkm2_completeness"]].max(axis=1) >= min_checkm
    ]
    print(f"Genome filtering: removing {len(metadata) - len(filtered)}/{len(metadata)} genomes")
    return filtered

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

def find_metaspecies(metadata: pd.DataFrame, min_size: int) -> dict[str, GTDBMetaSpecies]:
    """Find all metaspecies. Keeps only metaspecies with more than 'min_size' genomes"""
    species = metadata["gtdb_taxonomy"].str.split(";").str[-1].str.strip()
    metamap = defaultdict(lambda: defaultdict(int))

    for s in species:
        prefix = rm_suffix(s)
        metamap[prefix][s] += 1

    r =  {
        prefix: GTDBMetaSpecies(name=prefix, splits=dict(sc), genomes=defaultdict(list))
        for prefix, sc in metamap.items() if len(sc) > 1
    }

    for _, row in metadata.iterrows():
        sp = row["gtdb_taxonomy"].split(";")[-1].strip()
        meta_sp = rm_suffix(sp)
        if meta_sp in r:
            r[meta_sp].genomes[sp].append(row["genomes"])

    for k in list(r.keys()):
        if sum(len(x) for x in r[k].genomes.values()) < min_size:
            del r[k]

    return r

def filter_pangbank(metadata: pd.DataFrame, metaspecies: dict[str, GTDBMetaSpecies], min_size: int) -> tuple[pd.DataFrame, list[str]]:
    """Keep all non-split species with more than 'min_size' genomes, i.e. original pangbank filter.
    Keep all metaspecies with more than 'min_size' genomes, i.e. merge candidates.
    """
    species = metadata["gtdb_taxonomy"].str.split(";").str[-1].str.strip()
    species_counts = species.value_counts()

    ok_species = set(species_counts[species_counts >= min_size].index)
    ok_metaspecies = {name for name, ms in metaspecies.items() if ms.ok(min_size)}

    ok_species_no_splits = [s for s in ok_species if rm_suffix(s) not in metaspecies]

    print(f"#pangenomes (>= {min_size} genomes): {len(ok_species)}, {len(ok_species_no_splits)} without splits")
    print(f"#metaspecies (>= {min_size} genomes): {len(metaspecies)}")

    prefixes = species.map(rm_suffix)
    keep = species.isin(ok_species) | prefixes.isin(ok_metaspecies)

    return metadata[keep], ok_species_no_splits

def filter_input_genome(metadata: pd.DataFrame, genomes: set[str]) -> pd.DataFrame:
    filtered = metadata[metadata["genomes"].isin(genomes)]
    print(f"Genome input filtering: removing {len(metadata) - len(filtered)}/{len(metadata)} genomes")
    return filtered

def filter(metadata_path: Path,
           genomes: Path,
           min_checkm_repr: float,
           min_checkm: float,
           min_genome_count: int) -> tuple[pd.DataFrame, dict[str, GTDBMetaSpecies], list[str]]:

    s_genomes = set()
    with open(genomes) as f:
        for line in f:
            line = line.strip()
            if line:
                s_genomes.add(line.split("\t")[0].strip())

    metadata = parse_gtdb_metadata(metadata_path)
    metadata = filter_input_genome(metadata, s_genomes)
    metadata = filter_representative(metadata, min_checkm_repr)
    metadata = filter_genome(metadata, min_checkm)
    metaspecies = find_metaspecies(metadata, min_genome_count)
    metadata, ok_no_splits = filter_pangbank(metadata, metaspecies, min_genome_count)
    return metadata, metaspecies, ok_no_splits

def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Find split species using GTDB metadata files.",
        epilog="Example: python find_gtdb_splits.py --metadata-file <path> --ar-metadata-file <path>"
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
        "--used-genomes",
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
    path = args.metadata_file
    output = args.output_directory
    output.mkdir(parents=True, exist_ok=True)

    def process(path: Path):
        metadata, metaspecies, ok_no_splits = filter(
            path,
            args.used_genomes,
            args.representative_min_completeness,
            args.genome_min_completeness,
            args.min_genome_count,
        )
        out = output
        out.mkdir(exist_ok=True)
        metadata.to_csv(out / "meta.filtered.tsv", sep="\t")

        with open(out / "metaspecies.tsv", "w") as fo:
            for k, m in metaspecies.items():
                fo.write(f"{k}\t{';'.join(m.splits)}\n")

        pangenomes = {s : [] for s in ok_no_splits}
        for _, row in metadata.iterrows():
            s = row["gtdb_taxonomy"].split(";")[-1].strip()
            if s in pangenomes:
                pangenomes[s].append(row["genomes"])

        pangenome_out = out / "pangenomes"
        pangenome_out.mkdir(exist_ok=True)
        for name, genomes in pangenomes.items():
            fn = name.replace(" ", "_") + ".list"
            with open(pangenome_out / fn , "w") as fo:
                fo.writelines(s + "\n" for s in genomes)

        meta_out = out / "meta"
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

    process(path)

if __name__ == "__main__":
    sys.exit(main())
