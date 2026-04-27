#!/usr/bin/env Rscript

# classify_msi.R
#
# standalone mosaic classifier. reads all *.features.csv files from a
# directory, applies the two-threshold decision tree, writes a combined
# results csv.
#
# usage:
#   Rscript classify_msi.R \
#     --input_dir <dir> \
#     --pattern "*.features.csv" \
#     --threshold_high 0.0055 \
#     --threshold_low 0.0029 \
#     --output <results_csv>
#
# classifier logic (from mosaic_2.0_development.R lines 112-114):
#   peak_avg >= threshold_high                            -> MSI-H
#   threshold_low < peak_avg < threshold_high AND
#     defb_status == "unstable"                           -> MSI-H
#   otherwise                                             -> MSS
#
# samples with num_called == 0 emit msi_status = NA.

suppressPackageStartupMessages({
    library(data.table)
    library(optparse)
})

opts <- list(
    make_option("--input_dir",      type = "character", help = "directory with features csvs"),
    make_option("--pattern",        type = "character", default = "*.features.csv",
                help = "glob pattern for features files [default %default]"),
    make_option("--threshold_high", type = "double",    default = 0.0055,
                help = "peak_avg upper threshold [default %default]"),
    make_option("--threshold_low",  type = "double",    default = 0.0029,
                help = "peak_avg lower threshold [default %default]"),
    make_option("--output",         type = "character", help = "output results csv")
)
args <- parse_args(OptionParser(option_list = opts))

required <- c("input_dir", "output")
missing  <- required[vapply(required, function(k) is.null(args[[k]]), logical(1))]
if (length(missing)) stop("missing required args: ", paste(missing, collapse = ", "))

# glob -> regex (just ., *, ? — enough for typical patterns)
glob_to_regex <- function(g) {
    re <- gsub("\\.", "\\\\.", g)
    re <- gsub("\\*", ".*", re)
    re <- gsub("\\?", ".",   re)
    paste0("^", re, "$")
}
files <- list.files(args$input_dir, pattern = glob_to_regex(args$pattern), full.names = TRUE)
if (!length(files)) stop("no files matched '", args$pattern, "' in ", args$input_dir)

feats <- rbindlist(lapply(files, fread), use.names = TRUE, fill = TRUE)

classify_mosaic <- function(peak_avg, defb_status, num_called,
                            threshold_high = 0.0055,
                            threshold_low  = 0.0029) {
    if (!is.na(num_called) && num_called == 0) return(NA_character_)
    if (is.na(peak_avg))                       return(NA_character_)
    if (peak_avg >= threshold_high)            return("MSI-H")
    if (peak_avg >  threshold_low &&
        !is.na(defb_status) &&
        defb_status == "unstable")             return("MSI-H")
    "MSS"
}

feats[, msi_status := mapply(
    classify_mosaic,
    peak_avg, defb_status, num_called,
    MoreArgs = list(threshold_high = args$threshold_high,
                    threshold_low  = args$threshold_low)
)]

fwrite(feats, args$output, na = "NA")
