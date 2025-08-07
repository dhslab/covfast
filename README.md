covfast

A fast, C-based utility for calculating template-based coverage from CRAM/BAM files. covfast calculates depth of coverage statistics for specified genomic regions. It is designed to calculate template coverage, meaning that for a given read pair, any overlapping bases between the two reads are counted only once. This provides a measure of unique fragment coverage, which can be more accurate for certain applications than the standard alignment-based coverage. The program's logic and filtering are designed to closely match the output of samtools mpileup. 

Key Features:
- Template-Based Coverage: Avoids double-counting overlapping regions of read pairs. 
- High Performance: Written in C with htslib for efficient CRAM/BAM processing.
- Comparable to samtools: Mimics the filtering logic of samtools depth -s for comparable results.
- Flexible Filtering: Allows for custom setting of minimum mapping quality (MAPQ) and base quality.

Build Instructions:

Requirements: gcc make wget
Standard htslib dependencies: zlib, bzip2, lzma, curl, openssl.
Compiling: The included Makefile will automatically download and compile the required version of htslib before building the covfast executable.
Build the executable: make

This will create the covfast binary in the current directory.

Install to system (Optional):

To install the covfast executable to /usr/local/bin for system-wide access, run: sudo make install
Usage

Usage: ./covfast [OPTIONS] <bed_file> <cram_or_bam_file>

Positional Arguments:
  <bed_file>            A standard BED file defining regions of interest (0-based).
  <cram_or_bam_file>    An indexed CRAM or BAM file with read alignments.

Options:
  -r, --reference <file>    The reference genome FASTA file, required for CRAM.
  -o, --outfile <file>      Output file for statistics. [Default: stdout]
  -q, --minmapqual <int>    Minimum mapping quality. [Default: 1]
  -Q, --minbasequal <int>   Minimum base quality. [Default: 13]
  -c, --coverage_values <str> Comma-separated coverage thresholds for summary stats.
  -h, --help                Show this help message and exit.

Example:

covfast --reference /path/to/hg38.fa \
        --outfile coverage_results.tsv \
        /path/to/my_regions.bed \
        /path/to/sample.cram

Docker:

For a fully containerized and reproducible environment, you can use the provided Dockerfile. See the Docker README for detailed instructions.
