#! /usr/bin/env Rscript

# The Parker Lab
# theparkerlab.org
#
# Original author: Peter Orchard
# Modified by: Vivek Rai

suppressPackageStartupMessages(library(optparse))

option_list <- list(
  make_option(c("--peaks-dir"),
    action = "store", dest = "peaks_dir", type = "character",
    help = "[Required] Path to the directory of peak files ('*_peaks.broadPeak.noblacklist'"
  ),
  make_option(c("--number-samples-threshold"),
    action = "store", dest = "number_samples",
    default = 1, type = "numeric", help = "[Optional] Number of samples that must show overlap with a master peak"
  ),
  make_option(c("--fdr"),
    action = "store", type = "numeric", default = 0.05,
    help = "[Optional] FDR to filter initial peak lists by (default: 0.05)"
  ),
  make_option(c("--out"),
    action = "store", type = "character",
    help = "[Required] Name of the output file (a bed file of master peaks)"
  )
)

option_parser <- OptionParser(
  usage = "usage: Rscript %prog [options]",
  option_list = option_list, add_help_option = T
)
opts <- parse_args(option_parser)

suppressPackageStartupMessages({
  library("dplyr")
  library("tidyr")
  library("glue")
  library("GenomicRanges")
})

#
# Helper functions
#
parse_file_name <- function(f) {
  regex <- "^(.*)_peaks.broadPeak.noblacklist"
  lib <- gsub(regex, "\\1", basename(f))

  x <- c("library" = lib)
}


read_peak_file <- function(f) {
  print(glue("Reading {basename(f)}.."), stderr())

  tmp <- data.table::fread(f, header = F, stringsAsFactors = F, sep = "\t")
  if (ncol(tmp) != 9) {
    print("Pleae provide 9-column broadPeak file(s) as <FILE>.broadPeak.noblacklist.", stderr())
    exit(1)
  }
  colnames(tmp) <- c("chrom", "start", "end", "peak_name", "V5", "V6", "V7", "V8", "neg_log_10_q")

  tmp <- tmp %>%
    dplyr::filter(neg_log_10_q >= -log10(opts$fdr)) %>%
    dplyr::select(chrom:end, peak_name)

  tmp$library <- parse_file_name(f)["library"]
  tmp
}


#
# main
#
peaks_dir <- opts$peak_dir

PEAKS <- list.files(peaks_dir, pattern = "*.broadPeak.noblacklist", full.names = T)

results <- bind_rows(lapply(PEAKS, read_peak_file))

results <- results %>% filter(library != "EndoC_JI_2")


print(
  glue("Merging peaks from {length(unique(results$library))} libraries: {libraries}",
    libraries = paste(unique(results$library), collapse = ", ")
  ),
  stderr()
)

results.granges <- makeGRangesFromDataFrame(results, keep.extra.columns = T)
master_peaks <- GenomicRanges::reduce(results.granges)

print(glue("There are {n} master peaks", n = length(master_peaks)))

# Determine the overlap between master peaks and the peak calls from each library
overlaps <- findOverlaps(results.granges, master_peaks)
overlaps.master_peak <- master_peaks[subjectHits(overlaps)]
overlaps.library <- results$library[queryHits(overlaps)]

# visual sanity check
# Each row of the following objects should overlap:
# head(results[queryHits(overlaps),])
# head(overlaps.master_peak)

overlap_counts <- cbind(as.data.frame(overlaps.master_peak)[, 1:3], as.character(overlaps.library))
colnames(overlap_counts) <- c("chrom", "start", "end", "library")

o_master_peaks <- overlap_counts %>%
  unique() %>%
  dplyr::group_by(chrom, start, end) %>%
  dplyr::summarize(number_libraries_with_overlap = n()) %>%
  unique()

ggplot(o_master_peaks) + geom_bar(aes(number_libraries_with_overlap))
# ggplot(o_master_peaks) + geom_bar(aes(number_libraries_with_overlap, cumsum(..count..)))

o_master_peaks$chrom <- as.character(o_master_peaks$chrom)
o_master_peaks <- master_peaks[order(o_master_peaks$chrom, o_master_peaks$start), ]

print(glue("There are {n} master peaks", n = nrow(o_master_peaks)))

write.table(o_master_peaks,
  file = opts$out, append = F,
  quote = F, sep = "\t", row.names = F, col.names = F
)

print(glue("Overlap counts output to {out}", out = opts$out))
