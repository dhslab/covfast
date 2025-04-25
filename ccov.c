#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>
#include <ctype.h>
#include <math.h>

#include <htslib/sam.h>
#include <htslib/hts.h>
#include <htslib/faidx.h>
#include <htslib/bgzf.h>
#include <htslib/kstring.h>
#include <htslib/hts_defs.h>
#include <htslib/khash.h>

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
} args_t;

/*
 * parse_coverage_values():
 *   Split a comma-separated string of ints into an array.
 *   Returns the number of values parsed.
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
 *   Reads a line from a BGZF handle (gzipped or uncompressed).
 *   Uses htslib's kstring_t to handle variable-length lines.
 *   Returns 1 if a line was read, 0 on EOF or error.
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
 *   Given a single read, fill in rc->cov_mask for the bases that pass mapQ/baseQ
 *   and align within [rc->start, rc->end).  This function does a CIGAR walk to
 *   find matched reference positions. For each matched base, check base quality
 *   and set cov_mask[pos - start] = 1 if it passes.
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
            // Deletion or skipped region on the reference
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

/*
 * merge_read_cov():
 *   Combine rc2’s coverage bits into rc1.  Then free rc2.
 *   (Used if we want exactly one entry per read name in the hash.)
 */
static void merge_read_cov(read_cov_t *rc1, read_cov_t *rc2) {
    if (!rc1 || !rc2) return;
    if (rc1->start != rc2->start || rc1->end != rc2->end) {
        // Should never happen if region logic is consistent
        fprintf(stderr, "Warning: mismatch in read_cov start/end.\n");
        return;
    }
    for (int64_t i = 0; i < rc1->len; i++) {
        if (rc2->cov_mask[i]) {
            rc1->cov_mask[i] = 1;
        }
    }
    free(rc2->cov_mask);
    free(rc2);
}

/* --------------------------------------------------------------------------
   COVERAGE STATISTICS
   (Same as before, unchanged)
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
        fprintf(stderr, "ERROR: Memory allocation failed\n");
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
        fprintf(stderr, "ERROR: Memory allocation failed\n");
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
    // --------------------------------------------------------------------
    // Minimal argument parse (adapt as needed)
    if (argc < 9) {
        fprintf(stderr, "Usage: %s --bed_file <bed> --cram_file <cram> "
                        "--outfile <out> --reference <ref> "
                        "--minmapqual <int> --minbasequal <int> "
                        "[--coverage_values \"10,20,...\"]\n", argv[0]);
        return 1;
    }

    args_t args;
    memset(&args, 0, sizeof(args));

    // Default coverage values
    char cov_values_str[256] = "10,20,40,60,100,250,500,1000,1250,2500,4000";
    args.min_mapq = 1;
    args.min_baseq = 13;

    // Quick parse
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--bed_file") == 0) {
            args.bed_file = argv[++i];
        } else if (strcmp(argv[i], "--cram_file") == 0) {
            args.cram_file = argv[++i];
        } else if (strcmp(argv[i], "--outfile") == 0) {
            args.outfile = argv[++i];
        } else if (strcmp(argv[i], "--reference") == 0) {
            args.reference = argv[++i];
        } else if (strcmp(argv[i], "--minmapqual") == 0 || strcmp(argv[i], "-q") == 0) {
            args.min_mapq = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--minbasequal") == 0 || strcmp(argv[i], "-Q") == 0) {
            args.min_baseq = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--coverage_values") == 0 || strcmp(argv[i], "-c") == 0) {
            strncpy(cov_values_str, argv[++i], sizeof(cov_values_str)-1);
            cov_values_str[sizeof(cov_values_str)-1] = '\0';
        }
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

    // Read header
    bam_hdr_t *hdr = sam_hdr_read(in);
    if (!hdr) {
        fprintf(stderr, "ERROR: Cannot read header from %s\n", args.cram_file);
        sam_close(in);
        return 1;
    }

    // Set reference if given (important for CRAM or to enable ref-based checks)
    if (args.reference) {
        if (hts_set_fai_filename(in, args.reference) != 0) {
            fprintf(stderr, "ERROR: Failed to set reference %s\n", args.reference);
            bam_hdr_destroy(hdr);
            sam_close(in);
            return 1;
        }
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

    // --------------------------------------------------------------------
    // Process each BED line
    char linebuf[8192];
    while (read_bed_line(bedfp, linebuf, sizeof(linebuf))) {
        // Skip comments or empty lines
        if (linebuf[0] == '#' || strlen(linebuf) < 5) {
            continue;
        }

        // Expected format: chrom start end gene info
        char chrom[256], gene[256], info[256];
        int64_t start, end;
        if (sscanf(linebuf, "%255s %" SCNd64 " %" SCNd64 " %255s %255s",
                   chrom, &start, &end, gene, info) < 5) {
            // malformed or not enough fields
            continue;
        }
        int64_t region_len = end - start;
        if (region_len <= 0) {
            continue;
        }

        // Allocate final coverage array for [start, end)
        uint32_t *coverage = (uint32_t *)calloc(region_len, sizeof(uint32_t));
        if (!coverage) {
            fprintf(stderr, "ERROR: coverage calloc failed (region_len=%" PRId64 ")\n", region_len);
            continue;
        }

        // Setup a hash: read_name -> read_cov_t*
        khash_t(read2cov) *hmap = kh_init(read2cov);

        int tid = bam_name2id(hdr, chrom);
        if (tid < 0) {
            // Chrom not in header, skip
            free(coverage);
            kh_destroy(read2cov, hmap);
            continue;
        }

        // Create iterator for region
        hts_itr_t *iter = sam_itr_queryi(idx, tid, start, end);
        if (!iter) {
            free(coverage);
            kh_destroy(read2cov, hmap);
            continue;
        }

        // Read alignments in region
        bam1_t *b = bam_init1();
        while (sam_itr_next(in, iter, b) >= 0) {
            const char *read_name = bam_get_qname(b);
            if (!read_name) continue;

            // Insert or find existing entry in the hash
            int ret;
            khiter_t k = kh_put(read2cov, hmap, strdup(read_name), &ret);

            // If ret == 0, it means read_name was already in the hash.
            // If ret != 0, it's a new entry and we must allocate a new read_cov_t.
            if (ret != 0) {
                read_cov_t *rc = (read_cov_t*)calloc(1, sizeof(read_cov_t));
                rc->start   = start;
                rc->end     = end;
                rc->len     = region_len;
                rc->cov_mask= (char*)calloc(region_len, sizeof(char));
                kh_value(hmap, k) = rc;
            }

            // Now add coverage from this read
            read_cov_t *rc = kh_value(hmap, k);
            add_read_coverage(b, rc, args.min_mapq, args.min_baseq);
        }

        bam_destroy1(b);
        hts_itr_destroy(iter);

        // Now finalize coverage.  Each read_cov_t corresponds to one read name
        // (which may or may not have a mate). If the same read name appears
        // multiple times, we have merged coverage in the same read_cov_t struct.
        // We add '1' to the final coverage array for each position set in cov_mask.
        for (khiter_t k2 = kh_begin(hmap); k2 != kh_end(hmap); ++k2) {
            if (!kh_exist(hmap, k2)) continue;
            read_cov_t *rc = kh_value(hmap, k2);
            if (!rc) continue;
            for (int64_t i = 0; i < rc->len; i++) {
                if (rc->cov_mask[i]) {
                    coverage[i]++;
                }
            }
            free(rc->cov_mask);
            free(rc);

            // Because we used strdup on the read_name as the key,
            // we should free that key as well
            free((char*)kh_key(hmap, k2));
        }
        kh_destroy(read2cov, hmap);

        // Compute coverage stats and print
        compute_and_print_stats(chrom, start, end, gene, info,
                                coverage, region_len,
                                args.cov_values, args.n_cov_values,
                                outfp);

        free(coverage);
    }

    // Cleanup
    if (outfp && outfp != stdout) fclose(outfp);
    bgzf_close(bedfp);
    hts_idx_destroy(idx);
    bam_hdr_destroy(hdr);
    sam_close(in);

    return 0;
}
