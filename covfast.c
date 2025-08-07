/*
 * =====================================================================================
 *
 * Filename:  covfast.c
 *
 * Description:  Calculates template-based coverage for genomic regions.
 *
 * Version:  0.6
 * Created:  08/07/2024
 * Revision:  7
 * Compiler:  gcc
 *
 * Author:  dspencer
 *
 * =====================================================================================
 *
 * PROGRAM: covfast
 *
 * PURPOSE:
 * This program calculates depth of coverage statistics for specified genomic
 * regions from a CRAM or BAM file. It is designed to calculate "template"
 * coverage, meaning that for a given read pair (template), any overlapping
 * bases between the two reads are counted only once. This provides a measure
 * of unique fragment coverage, which can be more accurate for certain
 * applications than the standard alignment-based coverage calculated by tools
 * like 'samtools mpileup'.
 *
 * The program mimics the filtering and coverage logic of 'samtools mpileup',
 * including counting deletions and reference skips (N CIGAR operations) as
 * covered, to provide a comparable result.
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
 * - htslib (version 1.12 or later recommended).
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

#define COVFAST_VERSION "0.6"

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
   READ-BASED COVERAGE LOGIC
   We use a hash to collect partial coverage masks for each read name, ensuring
   that if two mates overlap, their union is counted exactly once.
   -------------------------------------------------------------------------- */

// A struct to store coverage for one read name (template) within [region_start, region_end).
typedef struct {
    int64_t start;         // region start
    int64_t end;           // region end
    int64_t len;           // end - start
    char   *cov_mask;      // cov_mask[i] = 1 if read covers (start + i), else 0
} read_cov_t;

// Create a hash from read-name -> read_cov_t*
KHASH_MAP_INIT_STR(read2cov, read_cov_t*)

/*
 * add_read_coverage():
 * Given a single read, fill in rc->cov_mask for the bases that pass mapQ/baseQ
 * and align within [rc->start, rc->end).  This function does a CIGAR walk to
 * find matched reference positions. For each matched base, check base quality
 * and set cov_mask[pos - start] = 1 if it passes.
 */
static void add_read_coverage(const bam1_t *b, read_cov_t *rc, int min_mapq, int min_baseq) {
    const bam1_core_t *c = &b->core;

    // If the read's mapping quality is below min_mapq, skip entirely.
    if (c->qual < min_mapq) return;

    // Get base qualities
    uint8_t *bq = bam_get_qual(b);

    // Start the CIGAR parsing
    uint32_t *cigar = bam_get_cigar(b);
    int64_t refpos = c->pos;  // reference position at start of this read alignment (0-based)
    int readpos = 0;          // read offset (0-based)

    for (int i = 0; i < c->n_cigar; i++) {
        int op  = bam_cigar_op(cigar[i]);
        int ol  = bam_cigar_oplen(cigar[i]);

        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
            // For each matched base, check baseQ, possibly set coverage
            for (int j = 0; j < ol; j++) {
                int64_t rpos = refpos + j;
                int rp = readpos + j;  // read offset for this base
                if (rp < 0 || rp >= c->l_qseq) continue; // sanity check

                if (bq[rp] >= min_baseq) {
                    if (rpos >= rc->start && rpos < rc->end) {
                        int64_t idx = rpos - rc->start;
                        rc->cov_mask[idx] = 1;
                    }
                }
            }
            refpos  += ol;
            readpos += ol;
        }
        else if (op == BAM_CDEL || op == BAM_CREF_SKIP) {
            // Deletion or reference skip. Count this as covered, similar to mpileup
            // Note: There is no base quality to check for a deleted base.
            for (int j = 0; j < ol; j++) {
                int64_t rpos = refpos + j;
                if (rpos >= rc->start && rpos < rc->end) {
                    int64_t idx = rpos - rc->start;
                    rc->cov_mask[idx] = 1;
                }
            }
            refpos += ol;
        }
        else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) {
            // Insertion or soft clip: consumes read bases, but not reference
            readpos += ol;
        }
        // BAM_CHARD_CLIP, BAM_CPAD, etc. do not consume either
        // so we do nothing for them.
    }
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

    // Sum, min, max
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

    // Sort for Q1, median, Q3
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

    // Compute percentage above thresholds
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

    // Print (tab-delimited)
    // Note: The internal 'start' is 0-based. We print start+1 to match 1-based BED convention.
    fprintf(outfp, "%s\t%" PRId64 "\t%" PRId64 "\t%s\t%s\t%" PRIu64 "\t%.1f\t%.1f\t%.1f\t%.1f\t%u\t%u",
            chrom, start + 1, end, gene, info,
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
    // --------------------------------------------------------------------
    // Argument parsing
    args_t args;
    memset(&args, 0, sizeof(args));

    // Default values
    char cov_values_str[256] = "10,20,40,60,100,250,500,1000,1250,2500,4000";
    args.min_mapq = 1;
    args.min_baseq = 13;
    args.threads = 1; // Default to 1 thread

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
            case 'h':
                print_usage(argv[0]);
                return 0;
            case 'v':
                printf("covfast version %s\n", COVFAST_VERSION);
                return 0;
            default:
                print_usage(argv[0]);
                return 1;
        }
    }

    // Check for positional arguments
    if (optind + 2 > argc) {
        fprintf(stderr, "Error: Missing positional arguments: <bed_file> and <cram_file>\n");
        print_usage(argv[0]);
        return 1;
    }
    args.bed_file = argv[optind];
    args.cram_file = argv[optind + 1];

    // Check for required reference argument
    if (!args.reference) {
        fprintf(stderr, "Error: --reference is a required argument.\n");
        print_usage(argv[0]);
        return 1;
    }

    // Parse coverage thresholds
    args.n_cov_values = parse_coverage_values(cov_values_str, args.cov_values, MAX_COV_VALUES);

    // --------------------------------------------------------------------
    // Open CRAM/BAM
    samFile *in = sam_open(args.cram_file, "r");
    if (!in) {
        fprintf(stderr, "ERROR: Cannot open CRAM/BAM file %s\n", args.cram_file);
        return 1;
    }

    // --- Set number of threads for CRAM/BAM decompression ---
    if (args.threads > 1) {
        hts_set_threads(in, args.threads);
    }

    // Read header
    bam_hdr_t *hdr = sam_hdr_read(in);
    if (!hdr) {
        fprintf(stderr, "ERROR: Cannot read header from %s\n", args.cram_file);
        sam_close(in);
        return 1;
    }

    // Set reference if given (important for CRAM or to enable ref-based checks)
    if (hts_set_fai_filename(in, args.reference) != 0) {
        fprintf(stderr, "ERROR: Failed to set reference %s\n", args.reference);
        bam_hdr_destroy(hdr);
        sam_close(in);
        return 1;
    }

    // Load index
    hts_idx_t *idx = sam_index_load(in, args.cram_file);
    if (!idx) {
        fprintf(stderr, "ERROR: Cannot load index for %s\n", args.cram_file);
        bam_hdr_destroy(hdr);
        sam_close(in);
        return 1;
    }

    // --------------------------------------------------------------------
    // Open BED file
    BGZF *bedfp = bgzf_open(args.bed_file, "r");
    if (!bedfp) {
        fprintf(stderr, "ERROR: Could not open BED file %s\n", args.bed_file);
        hts_idx_destroy(idx);
        bam_hdr_destroy(hdr);
        sam_close(in);
        return 1;
    }

    // Open outfile or use stdout
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

    // Print header line
    fprintf(outfp, "#chrom\tstart\tend\tgene\tinfo\ttotal_cvg\tmean_cvg\tQ1_cvg\tmedian_cvg\tQ3_cvg\tmin_cvg\tmax_cvg");
    for (int i = 0; i < args.n_cov_values; i++) {
        fprintf(outfp, "\tpct_above_%d", args.cov_values[i]);
    }
    fprintf(outfp, "\n");

    // --- Initialize reusable memory buffers ---
    uint32_t *coverage = NULL;
    int64_t coverage_buf_size = 0;
    khash_t(read2cov) *hmap = kh_init(read2cov);

    // --------------------------------------------------------------------
    // Process each BED line
    char linebuf[8192];
    while (read_bed_line(bedfp, linebuf, sizeof(linebuf))) {
        // Skip comments or empty lines
        if (linebuf[0] == '#' || strlen(linebuf) < 5) {
            continue;
        }

        // Expected format: chrom start end [gene] [info]
        char chrom[256], gene[256] = "N/A", info[256] = "N/A";
        int64_t start, end;
        // Use sscanf to be robust to missing optional fields
        int items_scanned = sscanf(linebuf, "%255s %" SCNd64 " %" SCNd64 " %255s %255s",
                                   chrom, &start, &end, gene, info);
        if (items_scanned < 3) {
            fprintf(stderr, "Warning: Skipping malformed BED line: %s\n", linebuf);
            continue;
        }

        // BED is 0-based, htslib iterator is 0-based. No conversion needed.
        int64_t region_len = end - start;
        if (region_len <= 0) {
            continue;
        }

        // --- Reuse or reallocate the coverage buffer ---
        if (region_len > coverage_buf_size) {
            coverage_buf_size = region_len;
            coverage = (uint32_t *)realloc(coverage, coverage_buf_size * sizeof(uint32_t));
            if (!coverage) {
                fprintf(stderr, "ERROR: coverage realloc failed\n");
                exit(1); // Critical memory failure
            }
        }
        // Always clear the buffer for the current region
        memset(coverage, 0, region_len * sizeof(uint32_t));


        int tid = bam_name2id(hdr, chrom);
        if (tid < 0) {
            fprintf(stderr, "Warning: Chromosome '%s' not found in CRAM header. Skipping region.\n", chrom);
            continue;
        }

        // Create iterator for region
        hts_itr_t *iter = sam_itr_queryi(idx, tid, start, end);
        if (!iter) {
            // This is not an error; it just means no reads are in the region.
        }

        // This bitmask contains all the flags we want to skip
        const uint32_t flags_to_skip = BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY;

        // Read alignments in region
        bam1_t *b = bam_init1();
        if (iter) {
            while (sam_itr_next(in, iter, b) >= 0) {
                if (b->core.flag & flags_to_skip) continue;

                const char *read_name = bam_get_qname(b);
                if (!read_name) continue;

                int ret;
                khiter_t k = kh_put(read2cov, hmap, strdup(read_name), &ret);

                if (ret != 0) {
                    read_cov_t *rc = (read_cov_t*)calloc(1, sizeof(read_cov_t));
                    rc->start   = start;
                    rc->end     = end;
                    rc->len     = region_len;
                    rc->cov_mask= (char*)calloc(region_len, sizeof(char));
                    kh_value(hmap, k) = rc;
                }

                read_cov_t *rc = kh_value(hmap, k);
                add_read_coverage(b, rc, args.min_mapq, args.min_baseq);
            }
        }

        bam_destroy1(b);
        if (iter) hts_itr_destroy(iter);

        // --- Finalize coverage (Pass 1) ---
        for (khiter_t k2 = kh_begin(hmap); k2 != kh_end(hmap); ++k2) {
            if (kh_exist(hmap, k2)) {
                read_cov_t *rc = kh_value(hmap, k2);
                for (int64_t i = 0; i < rc->len; i++) {
                    if (rc->cov_mask[i]) {
                        coverage[i]++;
                    }
                }
            }
        }

        // --- Free hash map contents (Pass 2) ---
        for (khiter_t k2 = kh_begin(hmap); k2 != kh_end(hmap); ++k2) {
            if (kh_exist(hmap, k2)) {
                read_cov_t *rc = kh_value(hmap, k2);
                free(rc->cov_mask);
                free(rc);
                free((char*)kh_key(hmap, k2));
            }
        }
        kh_clear(read2cov, hmap); // Reset the hash map for the next region
        
        // Compute coverage stats and print
        compute_and_print_stats(chrom, start, end, gene, info,
                                coverage, region_len,
                                args.cov_values, args.n_cov_values,
                                outfp);
    }

    // Cleanup
    free(coverage);
    kh_destroy(read2cov, hmap);
    if (outfp && outfp != stdout) fclose(outfp);
    bgzf_close(bedfp);
    hts_idx_destroy(idx);
    bam_hdr_destroy(hdr);
    sam_close(in);

    return 0;
}
