# scbirlab/nf-dock

![GitHub Workflow Status (with branch)](https://img.shields.io/github/actions/workflow/status/scbirlab/nf-dock/nf-test.yml)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.00.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](https://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

**scbirlab/nf-dock** is a Nextflow DSL2 pipeline for GPU-accelerated molecular docking of compound libraries against protein targets. It fetches AlphaFold protein structures, identifies binding pockets automatically, and docks all compounds against all targets in parallel using [GNINA](https://github.com/gnina/gnina), producing score matrices ready for downstream analysis.

**Table of contents**

- [Processing steps](#processing-steps)
- [Requirements](#requirements)
- [Quick start](#quick-start)
- [Inputs](#inputs)
- [Outputs](#outputs)
- [Issues, problems, suggestions](#issues-problems-suggestions)
- [Further help](#further-help)

## Processing steps

### Phase 1 — Prepare

1. **Fetch protein structures** (`FETCH_STRUCTURES`): Downloads AlphaFold2 structural models from the [EBI AlphaFold Database](https://alphafold.ebi.ac.uk/) for each target using its UniProt accession.

2. **Detect binding pocket** (`DETECT_POCKET`): Runs [fpocket](http://fpocket.sourceforge.net/) on each structure to identify the most druggable binding pocket. Extracts the pocket centroid and calculates a docking box with configurable padding. If no pocket is detected, falls back to blind docking over the whole protein.

3. **Prepare ligand library** (`PREPARE_LIGANDS`): Reads the input SDF file, adds hydrogens, generates 3D conformers (ETKDGv3 + MMFF optimisation) for any flat structures, and splits the library into chunks of 100 compounds for parallel docking.

### Phase 2 — Dock

4. **Dock** (`GNINA_DOCK`): Docks every ligand chunk against every protein pocket in parallel. Uses GNINA with Vina-based search and CNN-based rescoring to produce, per compound, the best-pose Vina score, CNNscore, and CNNaffinity. GPU acceleration is enabled by default.

5. **Aggregate scores** (`AGGREGATE_SCORES`): Collects all per-chunk score TSVs and merges them into two Parquet matrices — one in long format (one row per compound–protein pair) and one in wide format (compounds × proteins).

## Requirements

### Software

You need to have Nextflow and either Anaconda/Mamba, Singularity, or Docker installed on your system.

Nextflow ≥ 24.00.0 is required.

An NVIDIA GPU is recommended for docking. CPU-only mode is supported via `--gnina_gpu false`, but will be significantly slower.

#### First time using Nextflow?

If you're at the Crick or your shared cluster has it already installed, try:

```bash
module load Nextflow Singularity
```

Otherwise, if it's your first time using Nextflow on your system and you have Conda installed, you can install it using `conda`:

```bash
conda install -c bioconda nextflow
```

You may need to set the `NXF_HOME` environment variable. For example,

```bash
mkdir -p ~/.nextflow
export NXF_HOME=~/.nextflow
```

To make this a permanent change, you can do something like the following:

```bash
mkdir -p ~/.nextflow
echo "export NXF_HOME=~/.nextflow" >> ~/.bash_profile
source ~/.bash_profile
```

## Quick start

Make a [sample sheet (see below)](#sample-sheet), place your ligand SDF library in an `inputs/` folder, and optionally create a [`nextflow.config` file](#inputs) in the directory where you want the pipeline to run. Then run Nextflow.

```bash
nextflow run scbirlab/nf-dock \
    --sample_sheet inputs/sample-sheet.csv \
    --ligands_sdf library.sdf
```

Each time you run the pipeline after the first time, Nextflow will use a
locally-cached version which will not be automatically updated. If you want
to ensure that you're using the very latest version of the pipeline, use
the `-latest` flag.

```bash
nextflow run scbirlab/nf-dock -latest \
    --sample_sheet inputs/sample-sheet.csv \
    --ligands_sdf library.sdf
```

If you want to run a particular tagged version of the pipeline, such as `v0.0.1`, you can do so using

```bash
nextflow run scbirlab/nf-dock -r v0.0.1 \
    --sample_sheet inputs/sample-sheet.csv \
    --ligands_sdf library.sdf
```

For help, use `nextflow run scbirlab/nf-dock --help`.

The first time you run the pipeline for a project, the software dependencies
in `environment.yml` will be installed. This may take several minutes.

## Inputs

The following parameters are required:

- `sample_sheet`: Path to a CSV with information about the protein targets to be docked against.
- `ligands_sdf`: Filename of the compound library SDF file, relative to the `inputs` directory.

The following parameters have default values which can be overridden if necessary.

| Parameter | Default | Description |
| --------- | ------- | ----------- |
| `inputs` | `"inputs"` | Directory containing input files (must contain the `ligands_sdf` file). |
| `outputs` | `"output"` | Directory where pipeline outputs will be written. |
| `gnina_exhaustiveness` | `8` | GNINA search exhaustiveness. Higher values give more thorough sampling at the cost of speed. |
| `gnina_num_modes` | `9` | Number of binding poses generated per compound per target. Only the best pose (by CNNscore) is retained in the output. |
| `gnina_cnn` | `"fast"` | GNINA CNN scoring model. Options include `fast` (default) and `default`. |
| `box_padding` | `4.0` | Padding in Å added in each direction around the fpocket centroid to define the docking search box. |
| `gnina_cpus` | `4` | Number of CPU threads allocated per docking job. |
| `gnina_gpu` | `true` | Enable NVIDIA GPU acceleration for docking. Set to `false` for CPU-only mode. |

The parameters can be provided either in the `nextflow.config` file or on the `nextflow run` command.

Here is an example of the `nextflow.config` file:

```nextflow
params {
    sample_sheet = "inputs/sample-sheet.csv"
    ligands_sdf  = "library.sdf"

    gnina_exhaustiveness = 16
    gnina_gpu            = true
}
```

Alternatively, you can provide the parameters on the command line:

```bash
nextflow run scbirlab/nf-dock \
    --sample_sheet inputs/sample-sheet.csv \
    --ligands_sdf library.sdf \
    --gnina_exhaustiveness 16
```

### Sample sheet

The sample sheet is a CSV file listing the protein targets to dock against. One structure is fetched and docked per row.

The file must have a header with the column names below, and one line per target.

| Column | Required | Description |
| ------ | -------- | ----------- |
| `protein_id` | Yes | Unique identifier for this protein entry. Used to label all downstream outputs. |
| `pdb_id` | Yes | PDB accession code. Accepted but reserved for future use. |
| `chain` | No | Chain identifier to extract from the PDB structure. Accepted but reserved for future use. |
| `uniprot_id` | Yes | UniProt accession. Used to fetch the AlphaFold2 structural model from EBI. |
| `organism` | No | Source organism name. Carried through to outputs for annotation. |
| `gene_name` | No | Gene name. Carried through to outputs for annotation. |

Here is an example of the sample sheet:

| protein_id | pdb_id | chain | uniprot_id | organism | gene_name |
| ---------- | ------ | ----- | ---------- | -------- | --------- |
| ecoli_folA | 1RX2 | A | P0ABQ4 | Escherichia coli | folA |
| ecoli_murA | 1UAE | A | P0A749 | Escherichia coli | murA |
| ecoli_gyrB | 1KZN | A | P0AES6 | Escherichia coli | gyrB |
| ecoli_fabI | 1C14 | A | P0AEK4 | Escherichia coli | fabI |

### Ligand library

The ligand library must be provided as an SDF file. Both 2D and 3D input structures are accepted — the pipeline adds hydrogens and generates 3D conformers automatically using RDKit (ETKDGv3 + MMFF). The library is split internally into chunks of 100 compounds to maximise parallelism.

## Outputs

Outputs are saved in the directory specified by `--outputs` (`output` by default).

### `docking/`

| File | Description |
| ---- | ----------- |
| `score_matrix_long.parquet` | All docking scores in long format. One row per compound–protein pair. |
| `score_matrix_wide.parquet` | Score matrix in wide format. Compounds as rows, proteins as columns, `cnn_affinity` as values. |

The long-format Parquet file contains the following columns:

| Column | Description |
| ------ | ----------- |
| `compound_id` | Compound name (from SDF `_Name` field). |
| `protein_id` | Protein identifier (from sample sheet). |
| `uniprot_id` | UniProt accession (from sample sheet). |
| `vina_score` | AutoDock Vina score (kcal/mol) for the best pose. |
| `cnn_score` | GNINA CNN pose quality score (higher is better). |
| `cnn_affinity` | GNINA CNN predicted binding affinity (pKd units, higher is better). |

## Issues, problems, suggestions

If you run into problems not covered here, add to the
[issue tracker](https://www.github.com/scbirlab/nf-dock/issues).

## Further help

Here are the help pages of the software used by this pipeline.

- [GNINA](https://github.com/gnina/gnina)
- [fpocket](http://fpocket.sourceforge.net/)
- [RDKit](https://www.rdkit.org/docs/)
- [nextflow](https://www.nextflow.io/docs/latest/index.html)
