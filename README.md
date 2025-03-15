# RNA-seq Analysis Workflow

This Snakemake workflow automates the RNA-seq analysis pipeline outlined in the `pipeline.sh` script. It includes all steps from reference genome preparation to differential expression analysis.

## Features

- Automatic resource allocation based on system capabilities
- Modular design with separate rules for each step of the analysis
- Parallel execution of tasks where possible
- Reproducible analysis with clearly defined dependencies

## Requirements

- Conda environment as specified in `environment.yaml`
- Input FASTQ files in the `samples/` directory with naming format: `{sample}_1.fastq.gz` and `{sample}_2.fastq.gz`

## Directory Structure

```
.
├── config.yaml           # Configuration parameters
├── environment.yaml      # Conda environment specification
├── DEG.R                 # R script for differential expression analysis
├── Snakefile             # The workflow definition
├── samples/              # Directory containing raw FASTQ files
│   ├── sample1_1.fastq.gz
│   ├── sample1_2.fastq.gz
│   └── ...
└── output/               # Output directory (created by the workflow)
    ├── clean/            # Trimmed FASTQ files
    ├── mapped/           # Aligned BAM files
    ├── counts/           # Gene and transcript counts
    ├── logs/             # Log files and QC reports
    └── rdata/            # R data files with DEG results
```

## Usage

1. Place your paired-end FASTQ files in the `samples/` directory.

2. If needed, adjust parameters in the `config.yaml` file.

3. Run the workflow:

   ```bash
   snakemake --use-conda -j <cores>
   ```

   where `<cores>` is the number of CPU cores to use.

4. To generate a workflow visualization:

   ```bash
   snakemake --dag | dot -Tpng > dag.png
   ```

## Workflow Steps

1. **Create directories**: Sets up the necessary directory structure
2. **Download reference**: Downloads and prepares the human reference genome and annotation
3. **Build RSEM index**: Creates an index for RSEM/STAR alignment
4. **Trim reads**: Removes adapters and low-quality sequences using fastp
5. **Map reads**: Aligns reads to the reference genome using STAR
6. **Quantify transcripts**: Counts transcript expression using RSEM
7. **Run MultiQC**: Generates quality control reports
8. **Run DEG analysis**: Performs differential expression analysis using DESeq2 and ReactomePA

## Modifying the Workflow

- To add new steps or modify existing ones, edit the `Snakefile`
- To change parameters for tools, edit the `config.yaml` file
- To add new dependencies, update the `environment.yaml` file

## Troubleshooting

- Check the log files in `output/logs/` for errors
- If you encounter resource issues, adjust the `resources` section in `config.yaml`
- For conda environment issues, try recreating the environment: `conda env create -f environment.yaml --force`