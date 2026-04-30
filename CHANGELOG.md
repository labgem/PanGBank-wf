# labgem/pangbank: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## dev

### `Added`

- The original input genome list (before dereplication) is now saved to `genome_dereplication/<species>/original_input_genomes.tsv` for traceability.
- When `--genome_metadata` is provided, genomes are now sorted by completeness before cluster representative selection during dereplication. The completeness score is derived as `floor(max(checkm2_completeness, checkm_completeness))`. The full sort order is: `max_completeness` (desc), `L90` (asc), `L75` (asc), `L50` (asc), `auN` (desc).

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
