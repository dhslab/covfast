
## Overview

**`ccov`** is a C program that calculates per-base coverage statistics for genomic regions in a BED file using an indexed CRAM/BAM file.  
Unlike naive depth tools, it **merges overlapping mates** from paired-end reads, so each template base is counted only once.

Key features
------------

* Supports **CRAM, BAM, or SAM** (index required for CRAM/BAM).  
* Reads plain or **bgzipped (.bed.gz + .tbi) BED** files.  
* Filters by **mapping quality** (`-q`) and **base quality** (`-Q`).  
* Reports:
  - sum, mean, min, max, quartiles (Q1, median, Q3)  
  - % bases above custom **coverage thresholds** (defaults: 10,20,40,60,100,250,500,1000,1250,2500,4000)  
* Keeps only one BED region in memory — fast & memory-efficient.

---

## Requirements

* **htslib ≥ 1.18** (headers + library)  
* C11-compatible compiler (`gcc`, `clang`, …)  
* `make` (optional but convenient)

---

## Installation

```bash
git clone https://github.com/dhspence/ccov.git
cd ccov
gcc ccov.c -I htslib -L htslib/ -lhts -o ccov -Wl,-rpath="$(pwd)/htslib/" -lz -lpthread -lm

```

---

## Quick Start

```bash
./ccov \
  --bed_file     targets.bed.gz \
  --cram_file    sample.cram \
  --reference    GRCh38.fa \
  --outfile      sample.targets.cov.tsv \
  --minmapqual   20 \
  --minbasequal  13 \
  --coverage_values "10,30,100"
```

### Argument summary

| flag | required | description |
|------|----------|-------------|
| `--bed_file`        | ✔ | BED (plain or bgzipped) — 3/4/5 columns (`chrom start end [gene] [info]`) |
| `--cram_file`       | ✔ | Input CRAM/BAM/SAM **with index** (`.crai` / `.bai`) |
| `--reference`       |   | Reference FASTA (**mandatory for CRAM**) |
| `--outfile`         |   | Output TSV (default = stdout) |
| `--minmapqual, -q`  |   | Minimum mapping quality (default 1) |
| `--minbasequal, -Q` |   | Minimum base quality (default 13) |
| `--coverage_values, -c` |   | Comma-separated thresholds (≤ 50; default shown above) |

---

## Output columns

```
#chrom  start  end  gene  info  total_cvg  mean_cvg  Q1_cvg  median_cvg  Q3_cvg  min_cvg  max_cvg  pct_above_10 ...
```

* `chrom start end gene info` — copied from BED  
* `total_cvg` — sum of depths  
* `mean_cvg`, `Q1_cvg`, `median_cvg`, `Q3_cvg`, `min_cvg`, `max_cvg` — descriptive stats  
* `pct_above_X` — % positions ≥ X for each requested threshold

---

## Algorithm sketch

1. Iterate BED intervals.  
2. Build a **hash (read-name → bit-mask)** for the interval.  
3. For every alignment:
   * skip if `mapQ < minmapqual`; walk CIGAR.  
   * mark bases with `baseQ ≥ minbasequal`.  
   * mates share the same mask, so overlaps count once.  
4. Collapse masks into final depth array.  
5. Compute stats, write TSV line.  
6. Free memory, repeat.

---

## Example

```bash
bgzip -c regions.bed > regions.bed.gz
tabix -p bed regions.bed.gz

./ccov --bed_file regions.bed.gz \
          --cram_file NA12878.cram \
          --reference GRCh38.fa \
          --minmapqual 30 --minbasequal 20 \
          --outfile NA12878.chr7.cov.tsv
```

```bash
head NA12878.chr7.cov.tsv | column -t
```

---

## Contributing

Please open [issues](https://github.com/your-org/ccov/issues) or pull requests for bugs and feature requests.

---

## License

MIT License — see `LICENSE`.

---

## Citation

If `ccov` assists your work, cite this repository and **HTSlib**:

> Bonfield JK *et al.* HTSlib: C library for high-throughput sequencing data formats. *GigaScience* 2021.
```

