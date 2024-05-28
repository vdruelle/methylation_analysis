# methylation_analysis

## Introduction
This pipeline is used to analyze methylation patterns from `.bam` files generated using nanopore modified basecalling. It performs the following steps:
- Checks that the configuration file is provided.
- Checks that necessary input files exist relative to the config file.
- Aligns reads to the reference sequences using `dorado aligner`.
- Performs methylation pileup of alignments using `modkit`.
- Extracts methylation from regions/positions defined in BED files using `bedtools intersect`.

## How to run
### Prerequisite
Start by cloning this repo (on the cluster if aiming for cluster execution):
```
git clone https://github.com/vdruelle/methylation_analysis.git
```

Then create the conda environment for the pipeline using:
```
conda env create -f conda_env/methylation_analysis.yml
```

Additionally, the pipeline uses two specific tools which are expected to be available in the following locations:
- Dorado aligner, downloaded from https://github.com/nanoporetech/dorado and put in `softwares/dorado-0.7.0-linux-x64`
- Modkit, downloaded from https://github.com/nanoporetech/modkit and put in `softwares/modkit`

### Input format
The pipeline expects the following input directory structure:
```
input/
├── reads/
│ ├── reads1.bam
│ ├── reads2.bam
│ └── ...
├── references/
│ ├── ref1.fa
│ ├── ref2.fa
│ └── ...
└── regions/
├── ref1.bed
├── ref2.bed
└── ...
```
- `reads/`: Contains all the `.bam` files with reads.
- `references/`: Contains the reference sequences in fasta format, with the `.fa` extension.
- `regions/`: Contains the regions / positions on which to focus the analysis in `.bed` format. These files must have the same name as the fasta files of the references (except for the extension, which is `.bed` in this case). The format of this file has to be similar to:

```
CONTIG_NAME	start_pos1	end_pos1	.	0	DIRECTION
CONTIG_NAME	start_pos2	end_pos2	.	0	DIRECTION
```

The CONTIG_NAME is the name of the contig in the reference file (after `>` if the fasta file). The start and end positions define re regions to extract. If just focusing on single nucleotides, one can use pos 1000 and 1001 for example. The next fields, `.` and `0` need to stay like this. Finally DIRECTION defines in which strand direction to look for methylation patterns, can be `+`, `-` or `.` for both.

If looking for specific positions, it is suggested to look at the methylation pileup for the full genome to make sure the coordinates are right for the start_pos and end_pos.

### Config file
Create a configuration file named `run_config.yml` that corresponds to the data structure of your folder as described below:
```yaml
run_config:
  input: PATH_TO_INPUT
  output: PATH_TO_OUTPUT
  pileups:
    ref1_file:
      - "reads1"
      - "reads2"
    ref2_file:
      - "reads3"
      - "reads4"
```
You can save this file in the folder of your data for safekeeping. The `test_data` folder gives an example of the folder structure and its related config file.

### Running the pipeline
#### On the cluster (Scicore / SLURM) - recommended
To run the pipeline on scicore start by creating the input folder structure and config file as explained above. Then activate the environment using `conda activate methylation_analysis`. Then launch the pipeline using:
```bash
snakemake --profile cluster --configfile PATH_TO_CONFIG
```
Where `PATH_TO_CONFIG` is the path to your config file.

#### Local execution
To run the pipeline locally, activate the conda environment and use the command:
```bash
snakemake --configfile PATH_TO_CONFIG --cores NB_CORES
```
Replace `NB_CORES` by the number of cores you would like to allocate for the pipeline.

### Output
The pipeline generates the following output, all contained in the output folder as defined in your `run_config.yml`:

- Aligned `.bam` files and their indices in the `alignments` directory. These are organized per reference, as defined in the `run_config.yml`.
- Pileup files in the `pileup` directory. These are also organized per reference.
- Intersected BED files in the `focus_regions` directory. This is extracted from the full pileup for the positions of interest defined in the `regions` input folder.

The results, which in this case are the amount of methylation per genome position, is described in the files in the `pileup` and `focus_regions` folders. These files contain one line per base pair that was modified basecalled in the reference. Each of these lines has several columns and the content of each column is described in the table at the bottom of this documentation https://nanoporetech.github.io/modkit/intro_bedmethyl.html.

What is likely the most relevant for you is the position in column 2 and 3, the valid coverage in column 10 and the percentage of bases that were modified shown in column 11.