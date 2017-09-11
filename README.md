## Description

This tool is an attempt to do in silico PCRs from fastq or fasta files. The jellyfish version can only find perfect matches and takes about 3 times longer to get the same results as the bbduk version. It also only reports primer presence or absence. I left it for sentimental values!

**primer\_finder\_bbduk.sh** uses bbduk to find the primer sequences in the fastq/fasta file(s). If fastq file(s) are input, it pulls the matching reads out, which are then assembled with SPAdes. Then the primers are BLASTed on the assembly to try to see if both forward and reverse primers can be found in a single contig, thus a valid PCR product. PCR product length is also reported.

The longer the reads in the fastq file(s), the better the assembly and the less false negatives.

## Dependencies
1. bbduk (http://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/)
2. SPAdes (http://cab.spbu.ru/software/spades/)
3. BLAST+ (https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)

## Inputs

1. Primer pair list (fasta). Headers have to end with “_FWR” and “_REV”.
2. Raw reads (fastq) or assembly (fasta)
    * Single or multiple single-end reads
    * **One set of paired-end reads**
    * Single or multiple genomes (or assemblies)

## Usage

Typical command line using Illumina's paired-end fastq files:
````
bash primer_finder_bbduk.sh \
  -q -m \
  -n 1 \
  -p primers.fasta \
  -o output/ \
  mysample_R1.fastq.gz mysample_R2.fastq.gz

````

## Notes

The first version was based on Jellyfish. It worked well (just detects primer presence or absence), but when I tested bbduk, I was getting the exact same results in a third of the time. Besides speed (btw it takes between 1-2 min to process a pair of fastq files), the main advantage of using bbduk is that it allows up to 3 mismatches. I use 1 mismatch as default.
