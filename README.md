# pyMOI

## Introduction

This is a reimplementation of [estMOI](https://github.com/sammy-assefa/estMOI) in python. Still work in progress.

## Installation

```
pip install git+https://github.com/jodyphelan/pyMOI.git
```

## Usage

```
pyMOI --bam <input_bam> --vcf <input_vcf> --outfile <output_json_file>
```

## Help

```
usage: pyMOI [-h] --bam BAM --vcf VCF --outfile OUTFILE [--maxdist MAXDIST] [--min_count MIN_COUNT]

A simple command line tool for Moi

options:
  -h, --help            show this help message and exit
  --bam BAM             BAM (or cram) file
  --vcf VCF             VCF file
  --outfile OUTFILE     Name of output file
  --maxdist MAXDIST     Maximum distance between the first and last SNP
  --min_count MIN_COUNT
                        Minimum count of haplotype
```
