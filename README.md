## Description

This tool is an attempt to do _in silico_ PCRs from .fastq or .fasta files. 

The script ses bbduk to bait reads containing the primer sequences from the .fastq files. It 
performs a second round of baiting on the original .fastq files with the newly created baited
.fastq files in order to hopefully have sequence data to span the entire amplicon. Resulting
double-baited read files are assembled into contigs using SPAdes. 
The assemblies are BLASTed against the primer file to determine if both forward and reverse 
primers can be found in a single contig, thus a valid PCR product. PCR product length is also 
reported.

The longer the reads in the fastq file(s), the better the assembly and the fewer false negatives.

## External Dependencies
1. bbduk (http://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/)
2. SPAdes (http://cab.spbu.ru/software/spades/)
3. BLAST+ (https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
4. conda
5. pip


## Installation
Download the repository

`git clone https://github.com/adamkoziol/in_silico_PCR.git`

Enter the _in\_silico\_PCR_ directory
`cd in_silico_PCR`

Install the requirements using a combination of conda and pip

`while read requirement; do conda install --yes $requirement; done < requirements.txt`

`pip install OLCTools sipprverse biopython tables==3.4.2`

Add the in\_silico\_PCR folder to the $PATH

## Inputs

1. Primer pair list (fasta). Primer names have to end with “-F” or “-R”. Note: it is possible to have an integer 
following the direction: >vtx1a-F1 or >vtx1a-F are both acceptable
2. Raw reads (fastq) or assemblies (fasta)

## Usage

Typical command line usage:
````
primer_finder_bbduk.py ABSOLUTE PATH TO FOLDER IN WHICH REPORTS FOLDER TO BE CREATED -s ABSOLUTE PATH TO SEQUENCE FILES 
-p ABSOLUTE PATH AND NAME OF PRIMER FILE -m NUMBER OF MISMATCHES
````

## Test Dataset

I've provided six genomes in a mixture of FASTA and FASTQ formats to use to test the script on your system. 
The FASTA files are assemblies, while the FASTQ files are pre-baited files to reduce size. 
The report you create should match the one in the 'desired_outputs' folder (there may be small
differences when it comes to the order of genes).

NOTE: Please use absolute paths when running the program. If you don't there will probably be 
errors right at the beginning.

````
primer_finder_bbduk.py /path/to/in_silico_PCR/test_data -s /path/to/in_silico_PCR/test_data/sequences 
-p /path/to/in_silico_PCR/test_data/primers.txt
````

Note that in the provided ePCR.csv report file, in the results for _eae_ there are  four genome 
locations rather than the expected two. This is due to the fact that there are two separate 
primer sets for _eae_. This is fine for this output, but if you don't want results like this 
make sure that the names of the primers are unique.

Additionally, for 2014-SEQ-0121, the FASTQ and FASTA files yield different number of genes. The FASTQ sample has an
additional gene, _vtx2c_, this is likely due to issues with the assembly that do not always appear when
assembling with the double-baited reads (the de Bruijn graphs should be much simpler). 

## Options

````

usage: primer_finder_bbduk.py [-h] -s SEQUENCEPATH [-n CPUS] -p PRIMERFILE
                              [-m MISMATCHES]
                              path

Perform in silico PCR using bbduk and SPAdes

positional arguments:
  path                  Specify input directory

optional arguments:
  -h, --help            show this help message and exit
  -s SEQUENCEPATH, --sequencepath SEQUENCEPATH
                        Path of folder containing .fasta/.fastq(.gz) files to
                        process.
  -n CPUS, --cpus CPUS  Number of threads. Default is the number of cores in
                        the system
  -p PRIMERFILE, --primerfile PRIMERFILE
                        Absolute path and name of the primer file (in FASTA
                        format) to test. The file must haveevery primer on a
                        separate line AND -F/-R following the name e.g.
                        >primer1-F 
                        ATCGACTGACAC.... 
                        >primer1-R
                        ATCGATCGATCGATG.... 
                        >primer2-F 
                        .......
  -m MISMATCHES, --mismatches MISMATCHES
                        Number of mismatches allowed [0-3]. Default is 0

````
