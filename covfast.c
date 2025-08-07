/*
 * =====================================================================================
 *
 * Filename:  covfast.c
 *
 * Description:  Calculates template-based coverage for genomic regions.
 *
 * Version:  0.5
 * Created:  08/07/2024
 * Revision:  12
 * Compiler:  gcc
 *
 * Author:  David H. Spencer
 *
 * =====================================================================================
 *
 * PROGRAM: covfast
 *
 * PURPOSE:
 * This program calculates depth of coverage statistics for specified genomic
 * regions from a CRAM or BAM file. It is designed to calculate "template"
 * coverage by merging the coverage from read pairs, ensuring that overlapping
 * bases are counted only once. This is achieved using a memory-efficient
 * interval-based algorithm that is robust for very large regions.
 *
 * The program's filtering logic is designed to be comparable to standard tools
 * like 'samtools depth -s'.
 *
 * USAGE:
 * ./covfast [OPTIONS] <bed_file> <cram_or_bam_file>
 *
 * POSITIONAL ARGUMENTS:
 * bed_file            A standard BED file defining the regions of interest (0-based).
 * cram_or_bam_file    An indexed CRAM or BAM file with read alignments.
 *
 * OPTIONS:
 * --reference <ref.fa>    The reference genome FASTA file, required for CRAM.
 * Must be indexed with 'samtools faidx'. [Required]
 * --outfile <out.tsv>     Output file for statistics. [Default: stdout]
 * --threads, -t <int>     Number of threads for CRAM/BAM decompression. [Default: 1]
 * --minmapqual, -q <int>  Minimum mapping quality. [Default: 1]
 * --minbasequal, -Q <int> Minimum base quality. [Default: 13]
 * --coverage_values <str> Comma-separated coverage thresholds.
 * --help, -h              Show this help message and exit.
 * --version, -v           Show version number and exit.
 *
 * DEPENDENCIES:
 * - htslib (version 1.21 or later recommended).
 *
 * =====================================================================================
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>
#include <ctype.h>
#include <math.h>
#include <getopt.h> // For robust argument parsing

#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/faidx.h>
#include <htslib/bgzf.h>
#include <htslib/kstring.h>
#include <htslib/hts_defs.h>
#include <htslib/khash.h>

#define COVFAST_VERSION "0.5"

// Max coverage thresholds we can store
#define MAX_COV_VALUES 50

/* A structure to hold the parsed arguments */
typedef struct {
    char *bed_file;
    char *cram_file;
    char *outfile;
    char *reference;
    int min_mapq;
    int min_baseq;
    int n_cov_values;
    int cov_values[MAX_COV_VALUES];
    int threads; // Number of threads for I/O
} args_t;

void print_usage(char *prog_name) {
    fprintf(stderr, "\nUsage: %s [OPTIONS] <bed_file> <cram_or_bam_file>\n\n", prog_name);
    fprintf(stderr, "Calculates template-based coverage for genomic regions.\n\n");
    fprintf(stderr, "Positional Arguments:\n");
    fprintf(stderr, "  <bed_file>            A standard BED file defining regions of interest (0-based).\n");
    fprintf(stderr, "  <cram_or_bam_file>    An indexed CRAM or BAM file with read alignments.\n\n");
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -r, --reference <file>    The reference genome FASTA file, required for CRAM.\n");
    fprintf(stderr, "  -o, --outfile <file>      Output file for statistics. [Default: stdout]\n");
    fprintf(stderr, "  -t, --threads <int>       Number of threads for CRAM/BAM decompression. [Default: 1]\n");
    fprintf(stderr, "  -q, --minmapqual <int>    Minimum mapping quality. [Default: 1]\n");
    fprintf(stderr, "  -Q, --minbasequal <int>   Minimum base quality. [Default: 13]\n");
    fprintf(stderr, "  -c, --coverage_values <str> Comma-separated coverage thresholds for summary stats.\n");
    fprintf(stderr, "  -h, --help                Show this help message and exit.\n");
    fprintf(stderr, "  -v, --version             Show version number and exit.\n\n");
}


/*
 * parse_coverage_values():
 * Split a comma-separated string of ints into an array.
 * Returns the number of values parsed.
 */
int parse_coverage_values(const char *str, int *arr, int max_vals) {
    int count = 0;
    char *tmp = strdup(str);
    char *token = strtok(tmp, ",");
    while (token && count < max_vals) {
        arr[count++] = atoi(token);
        token = strtok(NULL, ",");
    }
    free(tmp);
    return count;
}

/*
 * read_bed_line():
 * Reads a line from a BGZF handle (gzipped or uncompressed).
 * Uses htslib's kstring_t to handle variable-length lines.
 * Returns 1 if a line was read, 0 on EOF or error.
 */
int read_bed_line(BGZF *fp, char *buf, size_t buflen) {
    static kstring_t str = {0, 0, NULL};

    int ret = bgzf_getline(fp, '\n', &str);
    if (ret < 0) {
        return 0;
    }
    strncpy(buf, str.s, buflen - 1);
    buf[buflen - 1] = '\0';
    return 1;
}

/*
 * Comparison function for qsort (ascending).
 */
int cmp_uint32_t(const void *a, const void *b) {
    uint32_t x = *((uint32_t *)a);
    uint32_t y = *((uint32_t *)b);
    if (x < y) return -1;
    else if (x > y) return 1;
    else return 0;
}

/* --------------------------------------------------------------------------
   READ-BASED COVERAGE LOGIC (INTERVAL-BASED)
   -------------------------------------------------------------------------- */

// A single covered interval [start, end)
typedef struct {
    int64_t start;
    int64_t end;
} cov_interval_t;

// A dynamic array of intervals for a single read template
typedef struct {
    size_t count;
    size_t capacity;
    cov_interval_t *intervals;
} read_cov_t;

// Create a hash from read-name -> read_cov_t*
KHASH_MAP_INIT_STR(read2cov, read_cov_t*)

// Comparison function for sorting intervals
int cmp_intervals(const void *a, const void *b) {
    const cov_interval_t *i1 = (const cov_interval_t *)a;
    const cov_interval_t *i2 = (const cov_interval_t *)b;
    if (i1->start < i2->start) return -1;
    if (i1->start > i2->start) return 1;
    return 0;
}

// Add a new interval to a read_cov_t struct
void add_interval(read_cov_t *rc, int64_t start, int64_t end) {
    if (start >= end) return; // Don't add empty intervals
    if (rc->count >= rc->capacity) {
        rc->capacity = rc->capacity == 0 ? 4 : rc->capacity * 2;
        rc->intervals = (cov_interval_t *)realloc(rc->intervals, rc->capacity * sizeof(cov_interval_t));
    }
    rc->intervals[rc->count].start = start;
    rc->intervals[rc->count].end = end;
    rc->count++;
}

// Merge sorted, overlapping intervals in place
void merge_intervals_inplace(read_cov_t *rc) {
    if (rc->count <= 1) return;

    qsort(rc->intervals, rc->count, sizeof(cov_interval_t), cmp_intervals);

    size_t merged_idx = 0;
    for (size_t i = 1; i < rc->count; i++) {
        // If current interval overlaps with or is adjacent to the last merged one, extend it
        if (rc->intervals[i].start <= rc->intervals[merged_idx].end) {
            if (rc->intervals[i].end > rc->intervals[merged_idx].end) {
                rc->intervals[merged_idx].end = rc->intervals[i].end;
            }
        } else {
            // No overlap, move to the next interval
            merged_idx++;
            rc->intervals[merged_idx] = rc->intervals[i];
        }
    }
    rc->count = merged_idx + 1;
}

/*
 * generate_read_intervals():
 * Creates a new read_cov_t struct containing the coverage intervals for a
 * single read, respecting the base quality filter.
 */
read_cov_t* generate_read_intervals(const bam1_t *b, int min_mapq, int min_baseq) {
    const bam1_core_t *c = &b->core;

    read_cov_t *rc = (read_cov_t *)calloc(1, sizeof(read_cov_t));
    uint8_t *bq = bam_get_qual(b);
    uint32_t *cigar = bam_get_cigar(b);
    int64_t refpos = c->pos;
    int readpos = 0;
    int64_t block_start = -1;

    for (int i = 0; i < c->n_cigar; i++) {
        int op = bam_cigar_op(cigar[i]);
        int ol = bam_cigar_oplen(cigar[i]);

        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
            for (int j = 0; j < ol; j++) {
                int rp = readpos + j;
                // Check if base quality passes
                if (rp >= 0 && rp < c->l_qseq && bq[rp] >= min_baseq) {
                    // If we are not in a block, start a new one
                    if (block_start == -1) {
                        block_start = refpos + j;
                    }
                } else {
                    // Base quality failed, end the current block if there is one
                    if (block_start != -1) {
                        add_interval(rc, block_start, refpos + j);
                        block_start = -1;
                    }
                }
            }
            refpos += ol;
            readpos += ol;
        } else if (op == BAM_CDEL || op == BAM_CREF_SKIP) {
            // End any existing block of matches
            if (block_start != -1) {
                add_interval(rc, block_start, refpos);
                block_start = -1;
            }
            // Deletions and ref skips do not contribute to coverage in a base-quality
            // sensitive model, as they have no base quality. Just advance refpos.
            refpos += ol;
        } else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) {
            // End any existing block, as these ops break contiguity on the reference
            if (block_start != -1) {
                add_interval(rc, block_start, refpos);
                block_start = -1;
            }
            readpos += ol;
        }
    }
    // Add the last block if the read ends with matches
    if (block_start != -1) {
        add_interval(rc, block_start, refpos);
    }
    
    // Merge any adjacent or overlapping intervals created during the process
    merge_intervals_inplace(rc);
    
    return rc;
}


/* --------------------------------------------------------------------------
   COVERAGE STATISTICS
   -------------------------------------------------------------------------- */
void compute_and_print_stats(const char *chrom, int64_t start, int64_t end,
                             const char *gene, const char *info,
                             uint32_t *coverage, int64_t region_len,
                             const int *cov_values, int n_cov_values,
                             FILE *outfp) {
    if (region_len <= 0) return;

    uint64_t sum_cov = 0;
    uint32_t min_cov = (uint32_t)(-1);
    uint32_t max_cov = 0;
    for (int64_t i = 0; i < region_len; i++) {
        uint32_t c = coverage[i];
        sum_cov += c;
        if (c < min_cov) min_cov = c;
        if (c > max_cov) max_cov = c;
    }
    double mean_cov = (double)sum_cov / (double)region_len;

    uint32_t *sorted = (uint32_t *)malloc(region_len * sizeof(uint32_t));
    if (!sorted) {
        fprintf(stderr, "ERROR: Memory allocation failed for sorted array\n");
        return;
    }
    memcpy(sorted, coverage, region_len * sizeof(uint32_t));
    qsort(sorted, region_len, sizeof(uint32_t), cmp_uint32_t);

    double q1=0.0, median=0.0, q3=0.0;
    if (region_len == 1) {
        q1 = median = q3 = (double)sorted[0];
    } else {
        int64_t idx_q1 = (int64_t)floor(0.25 * (region_len - 1));
        int64_t idx_q2 = (int64_t)floor(0.50 * (region_len - 1));
        int64_t idx_q3 = (int64_t)floor(0.75 * (region_len - 1));
        q1     = (double)sorted[idx_q1];
        median = (double)sorted[idx_q2];
        q3     = (double)sorted[idx_q3];
    }

    double *pct_above = (double *)calloc(n_cov_values, sizeof(double));
    if (!pct_above) {
        free(sorted);
        fprintf(stderr, "ERROR: Memory allocation failed for pct_above array\n");
        return;
    }

    for (int64_t i = 0; i < region_len; i++) {
        uint32_t c = coverage[i];
        for (int j = 0; j < n_cov_values; j++) {
            if (c > (uint32_t)cov_values[j]) {
                pct_above[j] += 1.0;
            }
        }
    }
    for (int j = 0; j < n_cov_values; j++) {
        pct_above[j] = 100.0 * pct_above[j] / (double)region_len;
    }

    fprintf(outfp, "%s\t%" PRId64 "\t%" PRId64 "\t%s\t%s\t%" PRIu64 "\t%.1f\t%.1f\t%.1f\t%.1f\t%u\t%u",
            chrom, start, end, gene, info,
            sum_cov, mean_cov, q1, median, q3, min_cov, max_cov);
    for (int j = 0; j < n_cov_values; j++) {
        fprintf(outfp, "\t%.1f", pct_above[j]);
    }
    fprintf(outfp, "\n");

    free(sorted);
    free(pct_above);
}

/* --------------------------------------------------------------------------
   MAIN
   -------------------------------------------------------------------------- */
int main(int argc, char *argv[]) {
    args_t args;
    memset(&args, 0, sizeof(args));

    char cov_values_str[256] = "10,20,40,60,100,250,500,1000,1250,2500,4000";
    args.min_mapq = 1;
    args.min_baseq = 13;
    args.threads = 1;

    static struct option long_options[] = {
        {"reference",       required_argument, 0, 'r'},
        {"outfile",         required_argument, 0, 'o'},
        {"threads",         required_argument, 0, 't'},
        {"minmapqual",      required_argument, 0, 'q'},
        {"minbasequal",     required_argument, 0, 'Q'},
        {"coverage_values", required_argument, 0, 'c'},
        {"help",            no_argument,       0, 'h'},
        {"version",         no_argument,       0, 'v'},
        {0, 0, 0, 0}
    };

    int c;
    while ((c = getopt_long(argc, argv, "r:o:t:q:Q:c:hv", long_options, NULL)) != -1) {
        switch (c) {
            case 'r': args.reference = optarg; break;
            case 'o': args.outfile = optarg; break;
            case 't': args.threads = atoi(optarg); break;
            case 'q': args.min_mapq = atoi(optarg); break;
            case 'Q': args.min_baseq = atoi(optarg); break;
            case 'c': strncpy(cov_values_str, optarg, sizeof(cov_values_str)-1); break;
            case 'h': print_usage(argv[0]); return 0;
            case 'v': printf("covfast version %s\n", COVFAST_VERSION); return 0;
            default: print_usage(argv[0]); return 1;
        }
    }

    if (optind + 2 > argc) {
        fprintf(stderr, "Error: Missing positional arguments: <bed_file> and <cram_file>\n");
        print_usage(argv[0]);
        return 1;
    }
    args.bed_file = argv[optind];
    args.cram_file = argv[optind + 1];

    if (!args.reference) {
        fprintf(stderr, "Error: --reference is a required argument.\n");
        print_usage(argv[0]);
        return 1;
    }

    args.n_cov_values = parse_coverage_values(cov_values_str, args.cov_values, MAX_COV_VALUES);

    samFile *in = sam_open(args.cram_file, "r");
    if (!in) {
        fprintf(stderr, "ERROR: Cannot open CRAM/BAM file %s\n", args.cram_file);
        return 1;
    }

    if (args.threads > 1) {
        hts_set_threads(in, args.threads);
    }

    bam_hdr_t *hdr = sam_hdr_read(in);
    if (!hdr) {
        fprintf(stderr, "ERROR: Cannot read header from %s\n", args.cram_file);
        sam_close(in);
        return 1;
    }

    if (hts_set_fai_filename(in, args.reference) != 0) {
        fprintf(stderr, "ERROR: Failed to set reference %s\n", args.reference);
        bam_hdr_destroy(hdr);
        sam_close(in);
        return 1;
    }

    hts_idx_t *idx = sam_index_load(in, args.cram_file);
    if (!idx) {
        fprintf(stderr, "ERROR: Cannot load index for %s\n", args.cram_file);
        bam_hdr_destroy(hdr);
        sam_close(in);
        return 1;
    }

    BGZF *bedfp = bgzf_open(args.bed_file, "r");
    if (!bedfp) {
        fprintf(stderr, "ERROR: Could not open BED file %s\n", args.bed_file);
        hts_idx_destroy(idx);
        bam_hdr_destroy(hdr);
        sam_close(in);
        return 1;
    }

    FILE *outfp = stdout;
    if (args.outfile && strcmp(args.outfile, "-") != 0) {
        outfp = fopen(args.outfile, "w");
        if (!outfp) {
            fprintf(stderr, "ERROR: Cannot open output file %s\n", args.outfile);
            bgzf_close(bedfp);
            hts_idx_destroy(idx);
            bam_hdr_destroy(hdr);
            sam_close(in);
            return 1;
        }
    }

    fprintf(outfp, "#chrom\tstart\tend\tgene\tinfo\ttotal_cvg\tmean_cvg\tQ1_cvg\tmedian_cvg\tQ3_cvg\tmin_cvg\tmax_cvg");
    for (int i = 0; i < args.n_cov_values; i++) {
        fprintf(outfp, "\tpct_above_%d", args.cov_values[i]);
    }
    fprintf(outfp, "\n");

    uint32_t *coverage = NULL;
    int64_t coverage_buf_size = 0;
    khash_t(read2cov) *hmap = kh_init(read2cov);

    char linebuf[8192];
    while (read_bed_line(bedfp, linebuf, sizeof(linebuf))) {
        if (linebuf[0] == '#' || strlen(linebuf) < 5) continue;

        char chrom[256], gene[256] = "N/A", info[256] = "N/A";
        int64_t start, end;
        int items_scanned = sscanf(linebuf, "%255s %" SCNd64 " %" SCNd64 " %255s %255s",
                                   chrom, &start, &end, gene, info);
        if (items_scanned < 3) {
            fprintf(stderr, "Warning: Skipping malformed BED line: %s\n", linebuf);
            continue;
        }

        int64_t region_len = end - start;
        if (region_len <= 0) continue;

        if (region_len > coverage_buf_size) {
            coverage_buf_size = region_len;
            coverage = (uint32_t *)realloc(coverage, coverage_buf_size * sizeof(uint32_t));
            if (!coverage) {
                fprintf(stderr, "ERROR: coverage realloc failed\n");
                exit(1);
            }
        }
        memset(coverage, 0, region_len * sizeof(uint32_t));

        int tid = bam_name2id(hdr, chrom);
        if (tid < 0) {
            fprintf(stderr, "Warning: Chromosome '%s' not found in CRAM header. Skipping region.\n", chrom);
            continue;
        }

        hts_itr_t *iter = sam_itr_queryi(idx, tid, start, end);
        const uint32_t flags_to_skip = BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY;

        bam1_t *b = bam_init1();
        if (iter) {
            while (sam_itr_next(in, iter, b) >= 0) {
                if (b->core.flag & flags_to_skip) continue;
                // skip if mapping quality is below threshold
                if (b->core.qual < args.min_mapq) continue;

                const char *read_name = bam_get_qname(b);
                if (!read_name) continue;

                read_cov_t *new_rc = generate_read_intervals(b, args.min_mapq, args.min_baseq);
                if (!new_rc) continue;

                int ret;
                // First, check if the key already exists without inserting.
                khiter_t k = kh_get(read2cov, hmap, read_name);

                if (k != kh_end(hmap)) { // Key already exists, merge.
                    read_cov_t *existing_rc = kh_value(hmap, k);
                    for(size_t i = 0; i < new_rc->count; ++i) {
                        add_interval(existing_rc, new_rc->intervals[i].start, new_rc->intervals[i].end);
                    }
                    merge_intervals_inplace(existing_rc);
                    // Free the new_rc object since its data has been merged.
                    free(new_rc->intervals);
                    free(new_rc);
                } else { // New key, so we insert it.
                    char *key_copy = strdup(read_name);
                    k = kh_put(read2cov, hmap, key_copy, &ret);
                    kh_value(hmap, k) = new_rc;
                }
            }
        }
        bam_destroy1(b);
        if (iter) hts_itr_destroy(iter);

        // Finalize coverage from intervals and free memory
        for (khiter_t k = kh_begin(hmap); k != kh_end(hmap); ++k) {
            if (kh_exist(hmap, k)) {
                read_cov_t *rc = kh_value(hmap, k);
                for (size_t i = 0; i < rc->count; ++i) {
                    for (int64_t pos = rc->intervals[i].start; pos < rc->intervals[i].end; ++pos) {
                        if (pos >= start && pos < end) {
                            coverage[pos - start]++;
                        }
                    }
                }
                free(rc->intervals);
                free(rc);
                free((char*)kh_key(hmap, k));
            }
        }
        kh_clear(read2cov, hmap);
        
        compute_and_print_stats(chrom, start, end, gene, info,
                                coverage, region_len,
                                args.cov_values, args.n_cov_values,
                                outfp);
    }

    free(coverage);
    kh_destroy(read2cov, hmap);
    if (outfp && outfp != stdout) fclose(outfp);
    bgzf_close(bedfp);
    hts_idx_destroy(idx);
    bam_hdr_destroy(hdr);
    sam_close(in);

    return 0;
}
