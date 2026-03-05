# scbirlab/nf-dock

![GitHub Workflow Status (with branch)](https://img.shields.io/github/actions/workflow/status/scbirlab/nf-dock/nf-test.yml)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.00.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](https://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

**scbirlab/nf-dock** is a Nextflow DSL2 pipeline for GPU-accelerated molecular docking of compound
libraries against protein targets. It fetches protein structures from RCSB PDB or AlphaFold,
identifies binding pockets automatically, and docks all compounds against all targets in parallel
using [GNINA](https://github.com/gnina/gnina), producing score matrices ready for downstream
analysis.

```
  sample_sheet.csv           ligands
  (protein targets)     (SDF or SMILES)
         │                    │
         ▼                    ▼
 FETCH_STRUCTURES        SplitLigands
  [PDB / AlphaFold]    [batches of 10]
         │                    │
         ▼                    ▼
  DETECT_POCKET        PREPARE_LIGANDS
  [fpocket top-3]    [add H, ETKDGv3,
         │             MMFF optimise]
         ▼                    │
   SplitPockets               │
  [pocket_00…02]              │
         │                    │
         └──────────┬─────────┘
                    ▼
              GNINA_DOCK
        [Vinardo + CNN rescore,
         sort by CNNaffinity]
                    │
                    ▼
      Extract_Gnina_scores
         [per-pose TSV]
                    │
                    ▼
        AGGREGATE_SCORES
   [score_matrix_long.tsv
    score_matrix_wide.tsv]
```

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

1. **`FETCH_STRUCTURES`**: Downloads a protein structure for each target.
   - If `pdb_id` is a plain accession (e.g. `1RX2`): fetches from RCSB PDB, extracts the specified
     chain, and runs pdbfixer to replace non-standard residues, add missing residues, and complete
     all heavy atoms.
   - If `pdb_id` begins with `AF2_` (e.g. `AF2_P0ABQ4`): fetches the corresponding AlphaFold2
     model from the EBI AlphaFold Database and completes heavy atoms with pdbfixer.

2. **`DETECT_POCKET`**: Runs [fpocket](http://fpocket.sourceforge.net/) on each structure and
   identifies up to the top 3 most druggable binding pockets. For each pocket it extracts the
   centroid coordinates and computes a bounding box (with configurable padding; minimum 20 Å per
   side). If no pocket is detected, falls back to a 30 Å blind-docking box centred on the origin.

3. **`SplitLigands`**: Reads the input compound library (SDF or SMILES) and splits it into batches
   of 10 compounds for parallel processing.

4. **`PREPARE_LIGANDS`**: For each batch, adds hydrogens, generates 3D conformers (ETKDGv3 + MMFF
   optimisation), and writes a prepared SDF. Molecules that fail conformer generation are skipped
   with a warning.

5. **`SplitPockets`**: Splits each protein's multi-pocket JSON into individual JSON files so that
   every pocket is docked independently.

### Phase 2 — Dock

6. **`GNINA_DOCK`**: Docks each prepared ligand batch against each protein pocket in parallel.
   Uses Vinardo scoring for pose generation and CNN rescoring to rerank poses by predicted
   binding affinity (`CNNaffinity`). Outputs are sorted by `CNNaffinity` and a docked SDF is
   written for every input. Runs on a dedicated GPU (`gpu_single` process label).

7. **`Extract_Gnina_scores`**: Parses the docked SDF output and writes a per-pose TSV containing
   ligand identity, receptor and pocket metadata, pose index, and all GNINA scores.

8. **`AGGREGATE_SCORES`**: Collects all per-batch TSV files and merges them into a long-format
   matrix (one row per pose) and a wide-format pivot table (best `cnn_affinity` per
   compound–protein pair).

## Requirements

### Software

You need to have Nextflow and either Anaconda/Mamba, Singularity, or Docker installed on your
system. Nextflow ≥ 24.00.0 is required.

An NVIDIA GPU is required for docking. In the `standard` profile, docking jobs run on a dedicated
`ga100` SLURM queue via the `gpu_single` process label; Singularity is launched with `--nv` so the
GPU is visible inside the container.

#### First time using Nextflow?

If you're at the Crick or your shared cluster has it already installed, try:

```bash
module load Nextflow Singularity
```

Otherwise, if it's your first time using Nextflow on your system and you have Conda installed, you
can install it using `conda`:

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

Make a [sample sheet (see below)](#sample-sheet), place your ligand library in the `inputs/`
directory, and optionally create a [`nextflow.config` file](#inputs) in the directory where you
want the pipeline to run. Then run Nextflow.

```bash
nextflow run scbirlab/nf-dock \
    --sample_sheet inputs/sample-sheet.csv \
    --ligands library.sdf
```

Each time you run the pipeline after the first time, Nextflow will use a locally-cached version
which will not be automatically updated. If you want to ensure that you're using the very latest
version of the pipeline, use the `-latest` flag.

```bash
nextflow run scbirlab/nf-dock -latest \
    --sample_sheet inputs/sample-sheet.csv \
    --ligands library.sdf
```

If you want to run a particular tagged version of the pipeline, such as `v0.0.1`, you can do so
using

```bash
nextflow run scbirlab/nf-dock -r v0.0.1 \
    --sample_sheet inputs/sample-sheet.csv \
    --ligands library.sdf
```

For help, use `nextflow run scbirlab/nf-dock --help`.

The first time you run the pipeline for a project, the software dependencies in `environment.yml`
will be installed. This may take several minutes.

## Inputs

The following parameters are required:

- `sample_sheet`: Path to a CSV with information about the protein targets to dock against.
- `ligands`: Filename of the compound library file within the `inputs` directory.

The following parameters have default values which can be overridden if necessary.

| Parameter | Default | Description |
| --------- | ------- | ----------- |
| `inputs` | `"inputs"` | Directory containing input files (must contain the `ligands` file). |
| `outputs` | `"output"` | Directory where all pipeline outputs will be written. |
| `batch_size` | `10` | Number of compounds per ligand batch. Larger batches reduce scheduling overhead; smaller batches increase parallelism. |
| `gnina_exhaustiveness` | `4` | GNINA search exhaustiveness. Higher values give more thorough sampling at the cost of speed. |
| `gnina_num_modes` | `3` | Number of binding poses generated per compound per pocket. All poses are retained in the output. |
| `gnina_cnn` | `"crossdock_default2018_KD_1"` | GNINA CNN model used for rescoring. |
| `box_padding` | `1.0` | Padding in Å added around the fpocket centroid in each direction. The box size is also clamped to a minimum of 20 Å per side. |
| `test` | `false` | Restrict the run to a small subset (5 proteins, 2 ligand batches, 3 docking jobs) to verify the pipeline runs end-to-end. |

The parameters can be provided either in the `nextflow.config` file or on the `nextflow run`
command.

Here is an example of the `nextflow.config` file:

```nextflow
params {
    sample_sheet = "inputs/sample-sheet.csv"
    ligands      = "library.sdf"

    gnina_exhaustiveness = 8
    gnina_cnn            = "crossdock_default2018_KD_1"
}
```

Alternatively, you can provide the parameters on the command line:

```bash
nextflow run scbirlab/nf-dock \
    --sample_sheet inputs/sample-sheet.csv \
    --ligands library.sdf \
    --gnina_exhaustiveness 8
```

### Sample sheet

The sample sheet is a CSV file listing the protein targets to dock against. One structure is
fetched and docked per row.

The file must have a header with the column names below, and one line per target.

| Column | Required | Description |
| ------ | -------- | ----------- |
| `protein_id` | Yes | Unique label for this target. Appears in all output filenames and score columns. |
| `pdb_id` | No | PDB accession (e.g. `1RX2`) to fetch from RCSB, or `AF2_<UniProt>` (e.g. `AF2_P0ABQ4`) to fetch the AlphaFold model instead. If blank, the AlphaFold model is fetched automatically using `uniprot_id`. |
| `chain` | No | Chain to extract from the PDB file (e.g. `A`). Ignored when fetching AlphaFold models. |
| `uniprot_id` | Yes | UniProt accession. Carried through to all score output columns. |
| `organism` | No | Source organism name. Carried through to outputs for annotation. |
| `gene_name` | No | Gene name. Carried through to outputs for annotation. |

Here is an example of the sample sheet:

| protein_id | pdb_id | chain | uniprot_id | organism | gene_name |
| ---------- | ------ | ----- | ---------- | -------- | --------- |
| ecoli_folA | 1RX2 | A | P0ABQ4 | Escherichia coli | folA |
| ecoli_murA | 1UAE | A | P0A749 | Escherichia coli | murA |
| ecoli_gyrB | 1KZN | A | P0AES6 | Escherichia coli | gyrB |
| ecoli_fabI | 1C14 | A | P0AEK4 | Escherichia coli | fabI |
| ecoli_ampC | 1FCO | A | P00811 | Escherichia coli | ampC |

### Ligand library

The compound library must be provided as either:

- **SDF** (2D or 3D) — molecules are read directly; 3D conformers are generated automatically.
- **SMILES** (`.smi`) — one SMILES string per line, with an optional name in the second column.

Hydrogens are added and 3D conformers are generated (ETKDGv3 + MMFF optimisation) automatically
for all input formats. The library is split into batches of 10 compounds to maximise parallelism.

## Outputs

Outputs are saved in the directory specified by `--outputs` (`output` by default).

### `receptors/`

Prepared PDB structure for each target, after chain extraction, residue fixing, and heavy-atom
completion. One file per protein, named `<protein_id>-<uniprot_id>-prepped.pdb`.

### `pockets/`

For each protein:

| File | Description |
| ---- | ----------- |
| `<protein_id>-<uniprot_id>-pocket.json` | JSON list of up to 3 pockets, each with centroid coordinates, box dimensions, and residue list. |
| `<protein_id>-<uniprot_id>-info.txt` | Raw fpocket summary text. |

### `ligands/`

Prepared 3D SDF files for each batch of compounds, named `<library>-mol.sdf`.

### `docking/`

Per-job docked pose SDF files and GNINA log files, plus the aggregated score matrices:

| File | Description |
| ---- | ----------- |
| `score_matrix_long.tsv` | All poses from all docking jobs. One row per pose. |
| `score_matrix_wide.tsv` | Best `cnn_affinity` per compound–protein pair, pivoted to a compounds × proteins matrix. |

Key columns in `score_matrix_long.tsv`:

| Column | Description |
| ------ | ----------- |
| `ligand_name` | Compound name from the input file (`_Name` SDF field or SMILES title). |
| `ligand_smiles` | SMILES string. |
| `ligand_zinc_id` | ZINC ID if present in the input file. |
| `receptor_uniprot_id` | UniProt accession of the target protein. |
| `receptor_pocket_id` | Pocket index (0, 1, or 2). |
| `receptor_pocket_center` | Pocket centroid coordinates [x, y, z]. |
| `receptor_pocket_size` | Docking box dimensions [x, y, z] in Å. |
| `blind_dock` | `True` if no pocket was detected and whole-protein docking was used instead. |
| `pose_id` | Sequential pose index within this docking job. |
| `library_chunk_id` | Index of the ligand batch that produced this pose. |
| `affinity` | Vinardo affinity score (kcal/mol) for this pose. |
| `cnn_score` | GNINA CNN pose quality score (higher is better). |
| `cnn_affinity` | GNINA CNN predicted binding affinity (pKd units, higher is better). |
| `cnn_affinity_var` | Variance of the CNN affinity prediction. |
| `cnn_vs` | GNINA CNN virtual screening score. |

## Issues, problems, suggestions

If you run into problems not covered here, add to the
[issue tracker](https://www.github.com/scbirlab/nf-dock/issues).

## Further help

Here are the help pages of the software used by this pipeline.

- [GNINA](https://github.com/gnina/gnina)
- [fpocket](http://fpocket.sourceforge.net/)
- [OpenBabel](https://openbabel.org/docs/)
- [pdbfixer](https://github.com/openmm/pdbfixer)
- [RDKit](https://www.rdkit.org/docs/)
- [nextflow](https://www.nextflow.io/docs/latest/index.html)
