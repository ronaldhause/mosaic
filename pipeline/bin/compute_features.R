#!/usr/bin/env Rscript

# compute_features.R
#
# per-sample feature extractor for the mosaic pipeline.
# reads an msings per-locus tsv for one sample and writes a single-row
# csv of sample-level features used by the mosaic classifier.
#
# usage:
#   Rscript compute_features.R \
#     --input <msings_tsv> \
#     --sample_name <name> \
#     --tumor_type <type> \
#     --num_loci 516876 \
#     --defb_locus "8:7679723-7679741" \
#     --output <features_csv>
#
# input formats accepted:
#   - full tcga format (20 cols incl. SAMPLE_NAME, LOCUS_COORDINATES,
#     PEAK_DIFFERENCE_VALUE, KS_VALUE, ...)
#   - simple 3-col format (sample, locus, peak_diff)
#
# output columns:
#   sample_name, tumor_type, peak_avg, peak_sd, num_unstable,
#   num_called, prop_unstable, defb_status

suppressPackageStartupMessages({
    library(data.table)
    library(optparse)
})

opts <- list(
    make_option("--input",       type = "character", help = "msings per-locus tsv"),
    make_option("--sample_name", type = "character", help = "sample name"),
    make_option("--tumor_type",  type = "character", default = NA, help = "tumor type label (optional)"),
    make_option("--num_loci",    type = "integer",   default = 516876L, help = "total loci in reference bed [default %default]"),
    make_option("--defb_locus",  type = "character", default = "8:7679723-7679741",
                help = "defb locus coordinates [default %default]"),
    make_option("--output",      type = "character", help = "output features csv")
)
args <- parse_args(OptionParser(option_list = opts))

required <- c("input", "sample_name", "output")
missing  <- required[vapply(required, function(k) is.null(args[[k]]), logical(1))]
if (length(missing)) stop("missing required args: ", paste(missing, collapse = ", "))

# load — data.table auto-handles tsv/csv; let it infer types then coerce
dat <- fread(args$input, sep = "\t", showProgress = FALSE)

# detect format. the 3-col format uses lowercase sample/locus/peak_diff;
# the tcga format uses SAMPLE_NAME / LOCUS_COORDINATES / PEAK_DIFFERENCE_VALUE.
cn <- colnames(dat)
if (all(c("LOCUS_COORDINATES", "PEAK_DIFFERENCE_VALUE") %in% cn)) {
    # canonical tcga format — nothing to do
} else if (all(c("locus", "peak_diff") %in% cn)) {
    setnames(dat, old = c("locus", "peak_diff"), new = c("LOCUS_COORDINATES", "PEAK_DIFFERENCE_VALUE"))
    if ("sample" %in% colnames(dat)) setnames(dat, "sample", "SAMPLE_NAME")
} else {
    stop("unrecognized input format. need either tcga-style (LOCUS_COORDINATES, PEAK_DIFFERENCE_VALUE) ",
         "or 3-col (sample, locus, peak_diff) columns. saw: ", paste(cn, collapse = ", "))
}

dat[, PEAK_DIFFERENCE_VALUE := as.numeric(PEAK_DIFFERENCE_VALUE)]

pdv          <- dat$PEAK_DIFFERENCE_VALUE
peak_avg     <- mean(pdv, na.rm = TRUE)
peak_sd      <- sd(pdv,   na.rm = TRUE)
num_unstable <- sum(pdv > 0, na.rm = TRUE)
num_na       <- sum(is.na(pdv))
num_called   <- args$num_loci - num_na
prop_unstable <- if (num_called > 0) num_unstable / num_called else NA_real_

# defb lookup — try the user-supplied coord, then with/without a leading "chr"
lookup_coords <- unique(c(
    args$defb_locus,
    sub("^chr", "", args$defb_locus),
    paste0("chr", sub("^chr", "", args$defb_locus))
))
defb_row <- dat[LOCUS_COORDINATES %in% lookup_coords]
defb_status <- if (nrow(defb_row) == 0 || is.na(defb_row$PEAK_DIFFERENCE_VALUE[1])) {
    "stable"
} else if (defb_row$PEAK_DIFFERENCE_VALUE[1] > 0) {
    "unstable"
} else {
    "stable"
}

out <- data.table(
    sample_name   = args$sample_name,
    tumor_type    = if (is.na(args$tumor_type)) NA_character_ else args$tumor_type,
    peak_avg      = peak_avg,
    peak_sd       = peak_sd,
    num_unstable  = num_unstable,
    num_called    = num_called,
    prop_unstable = prop_unstable,
    defb_status   = defb_status
)

fwrite(out, args$output)
