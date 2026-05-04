# labgem/pangbank: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## 0.1.0 - [2025-05-04]

### `Added`

- When `--genome_metadata` is provided, genomes are now sorted by completeness before cluster representative selection during dereplication. The completeness score is derived as `floor(max(checkm2_completeness, checkm_completeness))`. The full sort order is: `max_completeness` (desc), `L90` (asc), `L75` (asc), `L50` (asc), `auN` (desc).
- New `FILTER_GENOMES` process to filter input genomes based on CheckM/CheckM2 completeness scores. Activated with `--genome_quality_filtering`. Two thresholds are available: `--genome_min_completeness` (default: 70) for all input genomes and `--representative_min_completeness` (default: 85) for representative genomes to filter out species that have a poor quality representative genome.
- New `GTDB_SPLIT_SPECIES` subworkflow to detect and merge GTDB split species — species that GTDB splits but are genomically close enough to build a single pangenome. Activated with `--merge_gtdb_splits`. The ANI threshold for merging is controlled by `--gtdb_merge_ani_threshold` (default: 95).
- `PREPARE_PPANGGOLIN_INPUTS` (renamed from `PARSE_GENOMES_AND_TAXONOMY`) now handles merging of GTDB split species into metaspecies before writing per-species input files. Merged genomes are written to the metaspecies directory. A `pangenome_taxonomy.txt` file is written in each species output directory recording the pangenome-level taxonomy.
- Output files for genome preprocessing are now grouped under `genome_preprocessing/` with three sub-folders: `genome_filtering/`, `merge_gtdb_split/`, and `genome_dereplication/`.

### `Fixed`

- Fixed a non-deterministic publish race where the pre-dereplication `input_genomes.tsv.gz` could overwrite the post-dereplication one on `-resume`. A dedicated `PUBLISH_INPUT_GENOMES` step now writes the final file for all species from a single point in the DAG.

## 0.0.3 - [2025-04-01]

### `Added`

- Array specific profile inside custom config to improve overall cluster usage and increases scheduler friendliness
- Added support for TSV files containing genome-to-translation table mappings via the `--translation_tables` parameter. This allows to specify the correct translation table to use for a pangenome.

### `Changed`

- Update PPanGGOLiN version to 2.3.0

## 0.0.2 - [2025-07-08]

### `Added`

- You can now associate external metadata with genomes using the `--genome_metadata` parameter. These metadata are stored in the pangenome file.
- PPanGGOLiN process now produces tile plot to describe pangenome.
- Add genome in pangenome stats summary at pangenome level.
- Add proksee JSON map in the output.
- Compute genomes fluidity for each pangenome.

### `Fixed`

- Updated the storage of persistent families Mash sketch file to include only relative paths.

### `Dependencies`

### `Deprecated`

## 0.0.1 - [2025-01-14]

Initial release of labgem/pangbank, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- **Taxonomy and genomes processing** - Parse input genomes and the taxonomy to determine which species have enough genomes to build a pangenome
- **Dereplication** - Species with more genomes than a threshold are dereplicate using neighbor joining tree built from mash distances
- **Pangenome Construction** - Pangenomes are computed with PPanGGOLiN
- **Mash sketch of pangenome families** - Mash sektch of persistent families are built to be able to query easily a genome sequence and retrieve a matching pangenome.
