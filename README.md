# Nextflow pipeline template

![GitHub Workflow Status (with branch)](https://img.shields.io/github/actions/workflow/status/scbirlab/nf-dock/nf-test.yml)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.10.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](https://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

**scbirlab/nf-otemple** is a Nextflow pipeline template.

Use it to scaffold a new pipeline project.

**Table of contents**

- [Processing steps](#processing-steps)
- [Requirements](#requirements)
- [Quick start](#quick-start)
- [Inputs](#inputs)
- [Outputs](#outputs)
- [Issues, problems, suggestions](#issues-problems-suggestions)
- [Further help](#further-help)

## Processing steps


## Requirements

### Software

You need to have Nextflow and either Anaconda, Singularity, or Docker installed on your system.

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

Make a [sample sheet (see below)](#sample-sheet) and, optionally, 
a [`nextflow.config` file](#inputs) in the directory where you want the 
pipeline to run. Then run Nextflow.

```bash 
nextflow run scbirlab/nf-dock
```

Each time you run the pipeline after the first time, Nextflow will use a 
locally-cached version which will not be automatically updated. If you want 
to ensure that you're using the very latest version of the pipeline, use 
the `-latest` flag.

```bash 
nextflow run scbirlab/nf-dock -latest
```

If you want to run a particular tagged version of the pipeline, such as `v0.0.1`, you can do so using

```bash 
nextflow run scbirlab/nf-dock -r v0.0.1
```

For help, use `nextflow run scbirlab/nf-dock --help`.

The first time you run the pipeline for a project, the software dependencies 
in `environment.yml` will be installed. This may take several minutes.

## Inputs

The following parameters are required:

- `sample_sheet`: path to a CSV with information about the samples and FASTQ files to be processed

The following parameters have default values which can be overridden if necessary.

- `inputs = "inputs"` : The folder containing your inputs.
- `outputs = "outputs"` : The folder to containing the pipeline outputs.

The parameters can be provided either in the `nextflow.config` file or on the `nextflow run` command.

Here is an example of the `nextflow.config` file:

```nextflow
params {
    sample_sheet = "/path/to/sample-sheet.csv"
    inputs = "/path/to/inputs"
}
```

Alternatively, you can provide the parameters on the command line:

```bash
nextflow run scbirlab/nf-template \
    --sample_sheet /path/to/sample-sheet.csv \
    --inputs /path/to/inputs
``` 

### Sample sheet

The sample sheet is a CSV file providing information about which FASTQ files belong to which sample.

The file must have a header with the column names below, and one line per sample to be processed.

- `sample_id`: the unique name of the sample. 

Here is an example of the sample sheet:

| sample_id | 
| --------- |
| sample01  |
| sample02  |

## Outputs

Outputs are saved in the directory specified by `--outputs` (`outputs` by default). They are organised these directories:

- `multiqc`: HTML reports from the outputs of intermediate steps

## Issues, problems, suggestions

If you run into problems not covered here, add to the 
[issue tracker](https://www.github.com/scbirlab/nf-dock/issues).

## Further help

Here are the help pages of the software used by this pipeline.

- [multiqc](https://multiqc.info/)
- [nextflow](https://www.nextflow.io/docs/latest/index.html)
