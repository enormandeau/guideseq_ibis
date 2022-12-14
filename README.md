# GUIDEseq IBIS

GUIDE-seq analysis pipeline

Developed by [Eric Normandeau](https://github.com/enormandeau)

## TODO

* Modify `validate_project.sh`

- Start version system (only with tags on GitHub, no vx.y.z in text?)

## Description

GUIDE-seq IBIS is a pipeline for...

## Citation

GUIDE-seq IBIS is described in ...


## Use cases

- Research on cell cultures
- Agriculture plants
- ...

## Installation

To use GUIDE-seq IBIS, you will need a local copy of its repository. Different
releases can be [found here](https://github.com/enormandeau/guideseq_ibis/releases).
It is recommended to always use the latest release or even the developpment
version. You can either download an archive of the latest release at the above
link or get the latest commit (recommended) with the following git command:

```
git clone https://github.com/enormandeau/guideseq_ibis
```

### Dependencies

To run GUIDE-seq IBIS, you will also need to have the following programs installed
on your computer.

- GUIDE-seq IBIS will only work on GNU Linux or OSX
- bash 4+
- python 3.7+ (you can use miniconda3 to install python)
//- R 3+ (ubuntu/mint: `sudo apt-get install r-base-core`)
- java (ubuntu/mint: `sudo apt-get install default-jre`)
- [gnu parallel](https://www.gnu.org/software/parallel/)
- bwa 0.7.17-r1188+

### Preparation

- Install dependencies
- Download a copy of the GUIDE-seq IBIS repository (see **Installation** above)
- Get reference genome (ideally named `genome.fasta` in folder `03_genome`) and
  index it with `bwa index genome.fasta`
- Modify the parameters in `02_info/guideseq_ibis_config.sh` for your run
- Launch GUIDE-seq IBIS

## Overview of analyses

During the analyses, the following steps are performed:

- Filter and trim raw reads (`trimmomatic`)
- Extract UMI, dsODN, and other informations (Python script)
- Map reads to reference genome
- Remove PCR duplicates using UMI, first 8 nucleotides, and alignment position
- Find Double Strand Breaks (DSBs)
- Remove sites found in control
- Summarize results (Python script)

### Configuration file

Modify the parameters in `02_info/guideseq_ibis_config.sh` as needed.

### Launching the analysis

Launch the `guideseq_ibis` executable with the name of your configuration file as an
argument, like this:

```bash
./guideseq_ibis 02_info/guideseq_ibis_config.sh
```

### Running on the test dataset

If you want to test GUIDE-seq IBIS, jump straight to the `Test dataset` section
at the end of this file. Read through the README after to better understand the
program and it's outputs.

### Preparing samples

Copy your paired-end sample files in the `04_data` folder. You need one pair of
files per sample. The sequences in these files must contain the sequences of
the primer that you used during the PCR. Depending on the format in which you
received your sequences from the sequencing facility, you may have to proceed
to demultiplexing before you can use GUIDE-seq IBIS.

**IMPORTANT:** The file names must follow this format:

```
SampleID_*_R1_001.fastq.gz
SampleID_*_R2_001.fastq.gz
```

Notes: Each sample name, or SampleID, must contain no underscore (`_`) and be
followed by an underscore (`_`). The star (`*`) can be any string of text that
**does not contain space characters**. For example, you can use dashed (`-`) to
separate parts of your sample names, eg: `PopA-sample001_ANYTHING_R1_001.fastq.gz`.

### Results

Once the pipeline has finished running, all result files are found in the
`####_results` folder.

### Describe results by type

### Log files and parameters

For each run, three files are written in the `99_logfiles` folder. Each
contain a timestamp with the time of the run:

1. The exact config file that has been used
1. The full log of the run

## Test dataset

A test dataset is available as a [sister repository on
GitHub](https://github.com/enormandeau/guideseq_ibis_test_dataset). It is
composed of ####

If you have git and the GUIDE-seq IBIS dependencies installed, the following
commands will download the repository and the test data and put them in the
appropriate folder.

```bash
git clone https://github.com/enormandeau/guideseq_ibis
git clone https://github.com/enormandeau/guideseq_ibis_test_dataset
cp guideseq_ibis_test_dataset/04_data/* guideseq_ibis/04_data/
```

To run the analysis, move to the `guideseq_ibis` folder and launch:

```bash
cd guideseq_ibis
./guideseq_ibis 02_info/guideseq_ibis_config.sh
```

## License

CC share-alike

<a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/"><img alt="Creative Commons Licence" style="border-width:0" src="https://i.creativecommons.org/l/by-sa/4.0/88x31.png" /></a><br /><span xmlns:dct="http://purl.org/dc/terms/" property="dct:title">GUIDE-seq IBIS</span> by <span xmlns:cc="http://creativecommons.org/ns#" property="cc:attributionName">Eric Normandeau</span> is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/">Creative Commons Attribution-ShareAlike 4.0 International License</a>.<br />Based on a work at <a xmlns:dct="http://purl.org/dc/terms/" href="https://github.com/enormandeau/guideseq_ibis" rel="dct:source">https://github.com/enormandeau/guideseq_ibis</a>.
