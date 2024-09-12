[![GitHub Actions CI Status](https://github.com/PanGBank/workflows/nf-core%20CI/badge.svg)](https://github.com/PanGBank/actions?query=workflow%3A%22nf-core+CI%22)
[![GitHub Actions Linting Status](https://github.com/PanGBank/workflows/nf-core%20linting/badge.svg)](https://github.com/PanGBank/actions?query=workflow%3A%22nf-core+linting%22)

<!-- [![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX) -->

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

**PanGBank** is a bioinformatics pipeline that uses PPanGGOLiN to generate pangenomes from a list of input genomes and taxonomy. It then prepares files for the PanGBank API.

1. Parse input genomes and taxonomy, grouping them by species.
2. Run PPanGGOLiN on the genomes of each species.
3. Compile the pangenome information into a single file.
4. Create a mash sketch for all input genomes to easily identify the most appropriate pangenome for an input genome.
5. Concatenate the amino acid sequences of representative families into a single file for rapid retrieval of pangenomes containing a protein of interest.

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

PanGBank requires two input files:

1. **`--genomes <genome_file_list>`**
   A TSV file with two columns:

   - **Column 1:** `Genome_name` (unique name for each genome)
   - **Column 2:** Path to the corresponding genome file

2. **`--taxonomy <genome_taxonomy>`**
   A TSV file with two columns:
   - **Column 1:** `Genome_name` (must match the genome names in the `--genomes` file)
   - **Column 2:** Taxonomy, a list of taxon levels separated by a semicolon (`;`). The last taxon name is considered the species, and genomes will be grouped into species groups for pangenome analysis.

Now, you can run the pipeline using:

```bash
nextflow run labgem/pangbank \
   -profile <docker/singularity/.../institute> \
   --genomes <genome_file_list> \
   --taxonomy <genome_taxonomy>
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
> see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

## Credits

**PanGBank** was originally written in Snakemake and has been rewritten in Nextflow by Jean Mainguy.

<!-- We thank the following people for their extensive assistance in the development of this pipeline: -->

<!-- nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use PanGBank for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
