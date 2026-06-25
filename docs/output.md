# labgem/pangbank: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory after the pipeline has finished. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and build pangenomes from a input genomes using the following steps:

- [Genome quality filtering](#genome-quality-filtering) _(optional)_ - Filter input genomes based on CheckM/CheckM2 completeness scores
- [GTDB split species merging](#gtdb-split-species-merging) _(optional)_ - Detect and merge GTDB split species into metaspecies for pangenome construction
- [Taxonomy and genomes processing](#input-genomes-processing) - Parse input genomes and the taxonomy to determine which species have enough genomes to build a pangenome
- [Dereplication](#dereplication) - Species with more genomes than a threshold are dereplicated using a neighbor joining tree built from mash distances
- [Pangenome Construction](#pangenomes) - Pangenomes are computed with PPanGGOLiN
- [Mash sketch of pangenome families](#pangenomes) - Mash sketch of persistent families is built to allow querying a genome sequence and retrieving a matching pangenome
- [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

---

## Genome Preprocessing

All genome preprocessing outputs are written under `genome_preprocessing/`.

### Genome quality filtering

> Activated with `--genome_quality_filtering`

<details markdown="1">
<summary>Output files</summary>

- **`genome_preprocessing/genome_filtering/`**
  - **`genome_quality_filtering.info`**: Summary of the filtering step — number of genomes retained and removed, thresholds applied.
  - **`input_genomes.filtered.tsv`**: Filtered list of input genomes (passed to downstream steps).
  - **`genome_metadata.filtered.tsv`**: Filtered metadata file keeping only rows for retained genomes.

</details>

Genomes are filtered using completeness scores from CheckM or CheckM2 (taken from the `--genome_metadata` file). Two thresholds apply:

- `--genome_min_completeness` (default: 70) — minimum completeness for any input genome.
- `--representative_min_completeness` (default: 85) — minimum completeness for a genome to be eligible as a dereplication cluster representative.

### GTDB split species merging

> Activated with `--merge_gtdb_splits`

<details markdown="1">
<summary>Output files</summary>

- **`genome_preprocessing/merge_gtdb_split/`**
  - **`*.info`**: Per-species report from the split detection step, listing which species were identified as GTDB splits.
  - **`*.clusters`**: Per-species ANI clustering results used to determine which split species should be merged.
  - **`split_clusters.tsv`**: Concatenated cluster file passed to `PREPARE_PPANGGOLIN_INPUTS` to redirect merged genomes into a single metaspecies pangenome.
  - **`species_pair_summary.tsv`**: Per-species-pair ANI summary used to assess merge candidates. Each row represents a pair of species with columns: `mean_ani`, `mean_af`, `theoretical_pair_count`, `observed_count`, `passing_af_count`, `below_af_count` (pairs below `--gtdb_merge_af_threshold`, treated as ANI=0), and `missing_count` (pairs absent from the skani output, treated as ANI=0).
  </details>

GTDB sometimes splits a single biological species into several split species adding a letter suffix to the species name (ie : s\_\_Escherichia coli_F). This subworkflow clusters genomes from candidate split species by ANI (threshold: `--gtdb_merge_ani_threshold`, default: 95) and identifies groups that should be merged into a single pangenome (metaspecies). Genome pairs are only considered for merging if their alignment fraction meets `--gtdb_merge_af_threshold`; pairs below this threshold or absent from the skani output are treated as ANI=0.

### Input genomes processing

<details markdown="1">
<summary>Output files</summary>

- **`genome_preprocessing/`**
  - **`species_summary.tsv`**: Per-species count of genomes present in the taxonomy file and in the input dataset, and whether a pangenome will be built.
  - **`input_genomes/`**: Per-species directories containing the genome lists that will be used for dereplication or pangenome building.

</details>

### Dereplication

<details markdown="1">
<summary>Output files</summary>

- **`genome_preprocessing/genome_dereplication/<species>/`**
  - **`input_genomes.tsv.gz`**: Final list of dereplicated genomes passed to PPanGGOLiN.
  - **`original_input_genomes.tsv`**: Full input genome list before dereplication (for traceability).
  - **`mash_distance_matrix.phylip`**: Mash distance matrix used to build the neighbor joining tree.
  - **`nj_tree.nwk`**: Neighbor joining tree used for cluster representative selection.
  - **`clusters.tsv`**: Dereplication cluster assignments.

</details>

---

### Pangenomes

<details markdown="1">
<summary>Output files</summary>

- **`pangenome_summary.tsv`**: A summary file that compiles pangenome information for each species in a simple TSV format.

- **`pangenomes/`**
  - **`<species name>/`**: Directory containing the following files specific to the pangenome of a given species:
    - **`pangenome.h5`**: The main pangenome file generated by PPanGGOLiN, containing all data related to the pangenome.
    - **`pangenome_taxonomy.txt`**: The taxonomy string assigned to this pangenome. For standard species this is the full GTDB lineage. For metaspecies (GTDB split species merged together), this is the taxonomy of the 'metaspecie' without letter suffixes, since individual genomes may carry different species-level labels.
    - **`info.yaml`**: A YAML file with summary information about the pangenome.
    - **`input_genomes.tsv.gz`**: List of genomes and their path used to produce the pangenome.
    - **`genomes_md5sum.tsv.gz`**: A file listing the MD5 checksums for each input genome, used for version tracking.
    - **`genomes_statistics.tsv.gz`**: A compressed file summarizing the contents of each genome used in pangenome construction. Column descriptions are available in the PPanGGOLiN documentation: [Genome Statistics Table](https://ppanggolin.readthedocs.io/en/latest/user/PangenomeAnalyses/pangenomeAnalyses.html#genome-statistics-table).
    - **`persistent_nucleotide_families.fasta.gz`**: A compressed FASTA file containing the nucleotide sequences of persistent families.
    - **`all_protein_families.faa.gz`**: A compressed FASTA file containing the protein sequences of all families in the pangenome.
    - **`tile_plot.html`**: Tile plot of the pangenome. Visit PPanGGOLiN documentation for more detail : [U-shape plot documentation](https://ppanggolin.readthedocs.io/en/latest/user/PangenomeAnalyses/pangenomeAnalyses.html#tile-plot)
    - **`Ushaped_plot.html`**: Ushaped plot of families of the pangenome. Visit PPanGGOLiN documentation for more detail : [Tile plot documentation](https://ppanggolin.readthedocs.io/en/latest/user/PangenomeAnalyses/pangenomeAnalyses.html#u-shape-plot)
    - **`metadata/`**: A directory containing metadata associated with the pangenome. This metadata is also stored within the `pangenome.h5` file:
      - **`genomes_metadata_from_pangbank_wf_input.tsv.gz`**: A TSV file storing external genome metadata provided as input.
      - **`<genomes|contig>_metadata_from_annotation_file.tsv.gz`**: A TSV file containing metadata extracted from annotation files (GBFF or GFF) for genomes or contigs.
    - **`proksee/`**: A directory containing proksee JSON map for each genome of the pangenome.

</details>

[PPanGGOLiN](https://github.com/labgem/PPanGGOLiN) is a software suite used to create and manipulate prokaryotic pangenomes from a set of either genomic DNA sequences or provided genome annotations. For further reading and documentation see the [PPanGGOLiN documentation](https://ppanggolin.readthedocs.io/).

<!-- ### Dereplication -->

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reports generated by the pipeline: `pipeline_report.html`, `pipeline_report.txt` and `software_versions.yml`. The `pipeline_report*` files will only be present if the `--email` / `--email_on_fail` parameter's are used when running the pipeline.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
