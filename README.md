# GUIDE-Seq IBIS v0.1.0

## TODO

- Create test dataset (from the big dataset to come for the article, including annotation infos)
- Recommendation for `MIN_COVERAGE` parameter
- Modify `00_validate_project.sh`
- Bump to v1.0.0 before submission

## Description

GUIDE-Seq IBIS is used to validate CRISPR-Cas guide designs. In analyzes data
generated with guideseq experiments as modified by the [IBIS Genomic Analysis
Plateform](https://www.ibis.ulaval.ca/en/services-2/genomic-analysis-platform/).

Both the lab protocol and analysis pipeline were inspired by previous work
(https://github.com/tsailabSJ/guideseq)[https://github.com/tsailabSJ/guideseq],
although GUIDE-Seq IBIS has been completely re-written to accomodate the
protocol changes.

- Lab protocol: (GUIDE-Seq enables genome-wide profiling of off-target cleavage
by CRISPR-Cas nucleases)[https://pubmed.ncbi.nlm.nih.gov/25513782/]
- GUIDE-Seq IBIS repository on GitHub:
[https://github.com/enormandeau/guideseq_ibis](https://github.com/enormandeau/guideseq_ibis)


## Citation

GUIDE-Seq IBIS is described in ...

## Installation

To use GUIDE-Seq IBIS, you will need a local copy of its repository. Different
releases can be [found here](https://github.com/enormandeau/guideseq_ibis/tags).
It is recommended to always use the latest release or even the development
version. You can either download an archive of the latest release at the above
link or get the latest commit (recommended) with the following git command:

```
git clone https://github.com/enormandeau/guideseq_ibis
```

### Dependencies

To run GUIDE-Seq IBIS, you will also need to have the following programs installed
on your computer.

- GUIDE-Seq IBIS will only work on GNU Linux or OSX
- bash 4+
- python 3.7+ (you can use miniconda3 to install python)
- java (ubuntu/mint: `sudo apt-get install default-jre`)
- [gnu parallel](https://www.gnu.org/software/parallel/)
- bwa 0.7.17-r1188+
- OpenJDK (java)

### Preparation

- Install the dependencies
- Download a copy of the GUIDE-Seq IBIS repository
- Get reference genome (ideally named `genome.fasta` in folder `03_genome`)
- Index the genome, eg. with `bwa index genome.fasta`
- Add your paired-end GUIDE-Seq data to `04_data`
- Modify the parameters in `02_info/guideseq_ibis_config.sh` for your run
- Launch GUIDE-Seq IBIS

## Overview of analyses

During the analyses, the following steps are performed:

- Filter and trim raw reads (`trimmomatic`)
- Extract UMI, dsODN, and other informations (Python script)
- Map reads to reference genome (bwa)
- Remove PCR duplicates using UMI, first 8 nucleotides, and alignment position (Python script)
- Find Double Strand Breaks (DSBs, Python script)
- Summarize results (Python script)
- Optionally, annotate hits to show differences from expected guide

### Configuration file

Modify the parameters in `02_info/guideseq_ibis_config.sh` as needed.

### Launching the analysis

Launch the `guideseq_ibis` executable and pass the name of your configuration
file as an argument, like this:

```bash
./guideseq_ibis 02_info/guideseq_ibis_config.sh
```

### Running on the test dataset

If you want to test GUIDE-Seq IBIS, jump straight to the `Test dataset` section
at the end of this file. Read through the README after to better understand
GUIDE-Seq IBIS and its outputs.

### Preparing samples

Copy your paired-end sample files in the `04_data` folder. You need one pair of
files per sample. The sequences in these files must contain the sequences of
the primers that you used during the PCR.

**NOTE:** GUIDE-Seq IBIS does not demultiplex reads. Depending on the format in
which you received your sequences from the sequencing facility, you may have to
proceed to demultiplexing before you can use GUIDE-Seq IBIS.

**IMPORTANT:** The file names must follow this format:

```
SampleID_*_R1_001.fastq.gz
SampleID_*_R2_001.fastq.gz
```

Where `SampleID` is the name of the sample.

**NOTE**: Each sample name, or SampleID, must contain no underscore (`_`) and be
followed by an underscore (`_`). The star (`*`) can be any string of text that
**does not contain space characters**. For example, you can use dashes (`-`) to
separate parts of your sample names, eg: `GroupA-sample01_ANYTHING_R1_001.fastq.gz`.
In this case, the name of the sample is `GroupA-sample01`.

### Results

The results consist in:
1. A tsv file containing information about the enriched targets found for each
   sample

### Adding annotation

The output file can be annotated to show the differences between the sequences
of the off-target amplifications and the expected guide sequence:

```
./01_scripts/06_annotate_hits.py guideseq_ibis_report.tsv guide_info_per_sample.tsv \
    03_genome/genome.fasta guideseq_ibis_report_annotated.tsv
```

### Log files and parameters

For each run, two files are written in the `99_logfiles` folder. Each
contain a timestamp with the time of the run:

1. The exact config file that has been used
1. The full log of the run

## Test dataset

A test dataset is available as a [sister repository on
GitHub](https://github.com/enormandeau/guideseq_ibis_test_dataset). It is
composed of the first 250,000 reads from the samples used in the paper from
the `Citation` section near the top of this page.

If you have git and the GUIDE-Seq IBIS dependencies installed, the following
commands will download the repository and the test data and put them in the
appropriate folder.

1. Download GUIDE-Seq IBIS and the test dataset

```bash
git clone https://github.com/enormandeau/guideseq_ibis
git clone https://github.com/enormandeau/guideseq_ibis_test_dataset
cp guideseq_ibis_test_dataset/02_infos/* guideseq_ibis/02_infos
cp guideseq_ibis_test_dataset/04_data/* guideseq_ibis/04_data
cd guideseq_ibis
```

2. Download and index human genome

```
# Move to genome folder
cd 03_genome

# Get genome in archive
curl -OJX GET "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_000001405.40/download?include_annotation_type=GENOME_FASTA&filename=human_genome.zip" -H "Accept: application/zip"

# Extract it
unzip human_genome.zip

# Copy the assembly file
cp -l ncbi_dataset/data/GCF_000001405.40/GCF_000001405.40_GRCh38.p14_genomic.fna genome.fasta

# Index it
bwa index genome.fasta

# Move up to the main guideseq_ibis folder
cd ..
```

3. Run the analysis

```bash
# Run main analysis
./guideseq_ibis 02_info/guideseq_ibis_config.sh

# Annotation of the targets and off-targets
./01_scripts/06_annotate_hits.py guideseq_ibis_report.tsv guide_info_per_sample.tsv \
    03_genome/genome.fasta guideseq_ibis_report_annotated.tsv
```

## License

CC share-alike

<a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/"><img alt="Creative Commons Licence" style="border-width:0" src="https://i.creativecommons.org/l/by-sa/4.0/88x31.png" /></a><br /><span xmlns:dct="http://purl.org/dc/terms/" property="dct:title">GUIDE-Seq IBIS</span> by <span xmlns:cc="http://creativecommons.org/ns#" property="cc:attributionName">Eric Normandeau</span> is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/">Creative Commons Attribution-ShareAlike 4.0 International License</a>.
