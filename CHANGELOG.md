# labgem/pangbank: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## 0.0.2 - [date]

### `Added`

### `Fixed`

### `Dependencies`

### `Deprecated`

## 0.0.1 - [2025-01-14]

Initial release of labgem/pangbank, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- **Taxonomy and genomes processing** - Parse input genomes and the taxonomy to determine which species have enough genomes to build a pangenome
- **Dereplication** - Species with more genomes than a threshold are dereplicate using neighbor joining tree built from mash distances
- **Pangenome Construction** - Pangenomes are computed with PPanGGOLiN
- **Mash sketch of pangenome families** - Mash sektch of persistent families are built to be able to query easily a genome sequence and retrieve a matching pangenome.
