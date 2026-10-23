# LEAP

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A523.04.0-23aa8c.svg)](https://www.nextflow.io/)
[![Python](https://img.shields.io/badge/python-%E2%89%A53.7-3776AB?logo=python&logoColor=white)](https://www.python.org/)
[![Run with conda](https://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![Run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

LEAP is a Nextflow DSL2 pipeline for extending current gene annotations with long-read RNA-seq evidence. It identifies transcript-end extensions from long-read splice structures, validates 5' ends with CAGE data and 3' ends with polyA data, then applies safety checks to avoid creating unsupported overlaps with existing annotation.

The pipeline is organised around chromosome- and strand-aware processing so large BAM and annotation datasets can be run in parallel on HPC infrastructure.

## Contents

- [What LEAP Does](#what-leap-does)
- [Requirements](#requirements)
- [Inputs](#inputs)
- [Quick Start](#quick-start)
- [Parameters](#parameters)
- [Workflow](#workflow)
- [Outputs](#outputs)
- [Project Structure](#project-structure)
- [Configuration](#configuration)
- [Development Notes](#development-notes)
- [License](#license)

## What LEAP Does

- Converts long-read BAM alignments into BED, GenePred, and GTF-derived transcript structures.
- Extracts first-exon and last-exon unique splice sites from reads and current annotation.
- Uses CAGE and polyA BED files as independent evidence for transcript-end extensions.
- Retrieves current annotation from an Ensembl-style MySQL/MariaDB database.
- Filters partial transcripts and separates multi-exon from single-exon models.
- Detects natural overlaps in existing annotation so legitimate overlaps can be preserved.
- Selects the longest supported extension per terminal exon.
- Checks and fixes unsafe overlaps before generating extended transcript GFF files.

## Requirements

- Nextflow 23.04.0 or newer
- Linux or HPC execution environment
- Singularity/Apptainer or Conda
- Python 3 with `pandas`, `numpy`, and `matplotlib` available to helper scripts
- Perl and MySQL client tools for annotation database access
- Access to an annotation database compatible with `bin/db2gxf.pl`

Several modules declare Biocontainers images for tools such as `samtools`, `bedtools`, `bedToGenePred`, and `genePredToGtf`. The default config enables both Singularity and Conda:

```groovy
singularity.enabled = true
conda.enabled = true
```

## Inputs

LEAP expects:

- A coordinate-sorted long-read RNA-seq BAM file.
- A BAM index file (`.bai` or `.csi`).
- CAGE peaks in BED-like format for 5' validation.
- polyA sites in BED-like format for 3' validation.
- Database connection details for the current annotation source.
- An output directory for published results.

The CAGE and polyA files are split by chromosome and strand using the chromosomes found in the read BAM.

## Quick Start

Edit [params.config](/Users/rmerritt/my_app/LEAP/params.config) or provide parameters on the command line.

```groovy
params {
    output_directory = "/path/to/leap_results"

    reads = "/path/to/long_reads.bam"
    reads_index = "/path/to/long_reads.bam.bai"
    cage = "/path/to/cage.bed"
    polya = "/path/to/polya.bed"

    minimum_read_support = 2
    overlap_penalty = 10
    extend_on_current_annotation = true

    host = "database.host"
    port = 3306
    user = "username"
    pass = "password"
    dbname = "annotation_database"
}
```

Run locally with the default configuration:

```bash
nextflow run main.nf
```

Run on the EBI Slurm profile:

```bash
nextflow run main.nf -profile hpc_slurm
```

Resume a previous run:

```bash
nextflow run main.nf -profile hpc_slurm -resume
```

You can also override parameters at run time:

```bash
nextflow run main.nf \
  --reads /path/to/long_reads.bam \
  --reads_index /path/to/long_reads.bam.bai \
  --cage /path/to/cage.bed \
  --polya /path/to/polya.bed \
  --output_directory /path/to/leap_results \
  --minimum_read_support 2 \
  --overlap_penalty 10 \
  --extend_on_current_annotation true \
  --host database.host \
  --port 3306 \
  --user username \
  --pass password \
  --dbname annotation_database \
  -profile hpc_slurm
```

## Parameters

| Parameter | Required | Description |
| --- | --- | --- |
| `reads` | Yes | Long-read RNA-seq BAM file. |
| `reads_index` | Yes | Index for `reads`; used when extracting chromosomes and splitting BAMs. |
| `cage` | Yes | CAGE validation data. Used for first-exon/5' extension support. |
| `polya` | Yes | polyA validation data. Used for last-exon/3' extension support. |
| `output_directory` | Yes | Directory where selected pipeline outputs are published. |
| `minimum_read_support` | Yes | Minimum read support passed to read unique-splice-site extraction. |
| `overlap_penalty` | Yes | Alpha/penalty value passed to unsafe-overlap fixing. |
| `extend_on_current_annotation` | Yes | Controls whether candidate read extensions may be matched against current annotation. |
| `host` | Yes | Annotation database host. |
| `port` | Yes | Annotation database port. |
| `user` | Yes | Annotation database user. |
| `pass` | Yes | Annotation database password. |
| `dbname` | Yes | Annotation database name. |
| `sample_sheet` | No | Present in `params.config`, but not currently consumed by `main.nf`. |

`params.config` contains empty placeholders. Set values before running, either in that file or with `--parameter value` overrides.

## Workflow

The top-level workflow in [main.nf](/Users/rmerritt/my_app/LEAP/main.nf) runs four subworkflows.

### 1. `READ_PROCESSING`

Defined in [subworkflows/local/read_processing/main.nf](/Users/rmerritt/my_app/LEAP/subworkflows/local/read_processing/main.nf).

1. Extract chromosomes from the input BAM.
2. Split reads by chromosome.
3. Convert BAM to BED, BED to GenePred, and GenePred to GTF.
4. Split read-derived GTF files by strand.
5. Extract first and last exons.
6. Extract unique splice sites from read terminal exons.

### 2. `PROCESS_CAGE_AND_POLYA`

Defined in [subworkflows/local/process_cage_and_polya/main.nf](/Users/rmerritt/my_app/LEAP/subworkflows/local/process_cage_and_polya/main.nf).

1. Split CAGE data by chromosome and strand.
2. Split polyA data by chromosome and strand.
3. Emit keyed channels for downstream joins with read splice-site evidence.

### 3. `PROCESS_ANNOTATION`

Defined in [subworkflows/local/process_annotation/main.nf](/Users/rmerritt/my_app/LEAP/subworkflows/local/process_annotation/main.nf).

1. Retrieve current annotation per chromosome from the configured database.
2. Remove readthrough genes and OTT IDs during annotation retrieval.
3. Build transcript lists that should not be extended.
4. Filter partial transcripts.
5. Split multi-exon and single-exon transcripts.
6. Split current annotation by strand.
7. Extract first and last exons from current annotation.
8. Identify natural overlaps in current annotation.
9. Extract unique splice sites from current annotation.

### 4. `OVERLAP_READS`

Defined in [subworkflows/local/overlap_reads/main.nf](/Users/rmerritt/my_app/LEAP/subworkflows/local/overlap_reads/main.nf).

1. Overlap read-derived splice-site evidence with CAGE and polyA data.
2. Match validated candidate extensions against current annotation.
3. Select the longest possible extension.
4. Add gene IDs to selected extensions.
5. Filter internally overlapping features.
6. Intersect first-exon and last-exon extensions.
7. Check extension safety against natural overlaps.
8. Fix unsafe overlaps.
9. Build final extended transcript GFF files and per-chromosome statistics.

A more detailed workflow diagram is available in [BIOINFORMATICS_FLOWCHART.md](/Users/rmerritt/my_app/LEAP/BIOINFORMATICS_FLOWCHART.md).

## Outputs

Selected outputs are published under `params.output_directory`.

| Directory | Contents |
| --- | --- |
| `raw_annotation/` | Per-chromosome current annotation GFF3 files from the configured database. |
| `natural_overlaps/` | Existing annotation overlaps that are treated as natural/known overlaps. |
| `extension_safety/` | Safety-check outputs, including safe and unsafe extension BED files. |
| `unsafe_overlaps/` | Rejected unsafe overlaps, named like `<chromosome>_<strand>_rejected.bed`. |
| `fixed_extensions/` | First-exon and last-exon BED files after unsafe-overlap fixing. |
| `extended_transcripts/` | Final per-chromosome extended transcript GFF files. |
| `extended_transcripts_stats/` | Per-chromosome TSV statistics for the extended transcript output. |

Intermediate files remain in the Nextflow work directory unless a module publishes them.

## Project Structure

```text
LEAP/
├── main.nf
├── nextflow.config
├── params.config
├── conf/
│   └── ebi_slurm.config
├── bin/
│   ├── get_current_annotation.sh
│   ├── db2gxf.pl
│   └── *.py
├── assets/
│   └── *.py
├── modules/
│   └── local/
├── subworkflows/
│   └── local/
├── BIOINFORMATICS_FLOWCHART.md
├── LICENSE
└── README.md
```

## Configuration

The base [nextflow.config](/Users/rmerritt/my_app/LEAP/nextflow.config) includes [params.config](/Users/rmerritt/my_app/LEAP/params.config), enables Singularity and Conda, and defines the `hpc_slurm` profile.

The [conf/ebi_slurm.config](/Users/rmerritt/my_app/LEAP/conf/ebi_slurm.config) profile is configured for the EMBL-EBI Codon Slurm cluster. It sets the Slurm executor, queue behaviour, retry labels, and resource labels used by local modules.

## Development Notes

- The pipeline uses Nextflow DSL2 modules in `modules/local`.
- Reusable workflow units live in `subworkflows/local`.
- Helper scripts in `bin` are on the Nextflow process path and implement most of the annotation, overlap, and extension logic.
- Scripts in `assets` provide additional checking, filtering, and output-combination utilities.
- Keep module outputs keyed by chromosome, strand, and exon type when adding new steps; downstream joins depend on those keys.

## Author

Ryan Merritt

## License

See [LICENSE](/Users/rmerritt/my_app/LEAP/LICENSE).
