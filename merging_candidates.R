# merging_candidates.R
# Refactored by GitHub Copilot Chat Assistant for DanielOsmond
# Simplified and tidied: uses helper functions to reduce repetition

library(dplyr)
library(readr)
library(tidyr)
library(GenomicRanges)

# Configurable paths -------------------------------------------------------
root_dir <- "C:/Users/do271/OneDrive - University of Exeter/PhD/Genomics"
contigs_path <- file.path(root_dir, "contigs_names.txt")

# Helper functions ---------------------------------------------------------
read_contigs <- function(path) {
  contigs <- read.table(path, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  names(contigs) <- c("name1", "name2")
  contigs
}

read_and_parse_windows <- function(path, contigs, window_field = "window_id") {
  df <- read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  df <- df %>%
    separate(.data[[window_field]], into = c("chr", "pos"), sep = ":", remove = TRUE) %>%
    separate("pos", into = c("start", "end"), sep = "-", remove = TRUE)
  df <- merge(df, contigs, by.x = "chr", by.y = "name1", all.x = TRUE)
  if ("name2" %in% colnames(df)) df$chr <- df$name2
  df$chr <- suppressWarnings(as.numeric(df$chr))
  df$start <- as.numeric(df$start)
  df$end <- as.numeric(df$end)
  df
}

make_gr <- function(df, seqname_field = "chr") {
  makeGRangesFromDataFrame(
    df,
    keep.extra.columns = TRUE,
    seqinfo = NULL,
    seqnames.field = seqname_field,
    start.field = "start",
    end.field = "end",
    ignore.strand = TRUE,
    na.rm = TRUE
  )
}

# Overlap, optional filtering, and optional write --------------------------
overlap_and_maybe_write <- function(reference_gr, windows_df, out_path = NULL,
                                    seqname_field = "chr", filter_expr = NULL) {
  target_gr <- make_gr(windows_df, seqname_field)
  overlap_gr <- subsetByOverlaps(reference_gr, target_gr)
  out_df <- as.data.frame(overlap_gr)
  if (!is.null(filter_expr) && nzchar(filter_expr)) {
    out_df <- subset(out_df, eval(parse(text = filter_expr)))
  }
  if (!is.null(out_path)) write_tsv(out_df, out_path)
  out_df
}

# Read contigs --------------------------------------------------------------
contigs <- read_contigs(contigs_path)

# Read input window files --------------------------------------------------
Cu_eig1_99 <- read_and_parse_windows("./mafs/Cu_eig1_99%.txt", contigs)
Cu_eig1_99_fully <- read_and_parse_windows("./mafs/Cu_eig1_99%_fully_parallel.txt", contigs)

cymru_eig1_99 <- read_and_parse_windows("./mafs/cymru_eig1_99%.txt", contigs)
cymru_eig1_99_fully <- read_and_parse_windows("./mafs/cymru_eig1_99%_fully_parallel.txt", contigs)

NE_eig1_99 <- read_and_parse_windows("./mafs/NE_eig1_99%.txt", contigs)
NE_eig1_99_fully <- read_and_parse_windows("./mafs/NE_eig1_99%_fully_parallel.txt", contigs)

Pb_eig1_99 <- read_and_parse_windows("./mafs/Pb_eig1_99%.txt", contigs)
Pb_eig1_99_fully <- read_and_parse_windows("./mafs/Pb_eig1_99%_fully_parallel.txt", contigs)

UK_eig1_99 <- read_and_parse_windows("./mafs/UK_eig1_99%.txt", contigs)
UK_eig1_99_fully <- read_and_parse_windows("./mafs/UK_eig1_99%_fully_parallel.txt", contigs)

# Read picmin table and create ranges --------------------------------------
picmin05 <- read.table("./PicMin/picmin_fdr05_contribs.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
picmin05$start <- picmin05$mid - 25000
picmin05$end   <- picmin05$mid + 25000
picmin05 <- merge(picmin05, contigs, by.x = "chr", by.y = "name1", all.x = TRUE)
if ("name2" %in% colnames(picmin05)) picmin05$chr <- picmin05$name2
picmin05$chr <- suppressWarnings(as.numeric(picmin05$chr))
picmin05_gr <- make_gr(picmin05, seqname_field = "chr")

# PicMin overlaps ----------------------------------------------------------
# UK
overlap_and_maybe_write(picmin05_gr, UK_eig1_99, "UK.AF_picmin05_overlap.tsv")
overlap_and_maybe_write(picmin05_gr, UK_eig1_99_fully, "UK.AF.fully_picmin05_overlap.tsv")

# NE
overlap_and_maybe_write(picmin05_gr, NE_eig1_99, "NE.AF_picmin05_overlap.tsv",
                        filter_expr = "Tyne1_p <= 0.05 & Wear_p <= 0.05")
overlap_and_maybe_write(picmin05_gr, NE_eig1_99_fully, "NE.AF.fully_picmin05_overlap.tsv",
                        filter_expr = "Tyne1_p <= 0.05 & Wear_p <= 0.05")

# Cu
overlap_and_maybe_write(picmin05_gr, Cu_eig1_99, "Cu.AF_picmin05_overlap.tsv",
                        filter_expr = "Kernow_p <= 0.05 & Anglesey_p <= 0.05 & Ireland2_p <= 0.05")
overlap_and_maybe_write(picmin05_gr, Cu_eig1_99_fully, "Cu.AF.fully_picmin05_overlap.tsv",
                        filter_expr = "Kernow_p <= 0.05 & Anglesey_p <= 0.05 & Ireland2_p <= 0.05")

# Pb
overlap_and_maybe_write(picmin05_gr, Pb_eig1_99, "Pb.AF_picmin05_overlap.tsv",
                        filter_expr = "Nwales1_p <= 0.05 & Tyne1_p <= 0.05 & Wear_p <= 0.05 & MidWales_p <= 0.05")
overlap_and_maybe_write(picmin05_gr, Pb_eig1_99_fully, "Pb.AF.fully_picmin05_overlap.tsv",
                        filter_expr = "Nwales1_p <= 0.05 & Tyne1_p <= 0.05 & Wear_p <= 0.05 & MidWales_p <= 0.05")

# Cymru
overlap_and_maybe_write(picmin05_gr, cymru_eig1_99, "cymru.AF_picmin05_overlap.tsv",
                        filter_expr = "Nwales1_p <= 0.05 & Anglesey_p <= 0.05 & MidWales_p <= 0.05")
overlap_and_maybe_write(picmin05_gr, cymru_eig1_99_fully, "cymru.AF.fully_picmin05_overlap.tsv",
                        filter_expr = "Nwales1_p <= 0.05 & Anglesey_p <= 0.05 & MidWales_p <= 0.05")

# PBS overlaps --------------------------------------------------------------
pbs <- read_tsv("./sfs/pbs/gene_annotated_pbs.txt")

pbs_gr <- makeGRangesFromDataFrame(
  pbs,
  keep.extra.columns = TRUE,
  seqinfo = NULL,
  seqnames.field = "chromosome_name",
  start.field = "start",
  end.field = "end",
  ignore.strand = TRUE,
  na.rm = TRUE
)

# UK with PBS
overlap_and_maybe_write(pbs_gr, UK_eig1_99_fully, "UK.AF.fully_pbs_overlap.tsv")

# Cu with PBS
overlap_and_maybe_write(pbs_gr, Cu_eig1_99, "Cu.AF_pbs_overlap.tsv",
                        filter_expr = "AVB_pval <= 0.05 | HAY_pval <= 0.05 | AGA_pval <= 0.05")
overlap_and_maybe_write(pbs_gr, Cu_eig1_99_fully, "Cu.AF_fully_pbs_overlap.tsv",
                        filter_expr = "AVB_pval <= 0.05 | HAY_pval <= 0.05 | AGA_pval <= 0.05")

# Pb with PBS
overlap_and_maybe_write(pbs_gr, Pb_eig1_99, "Pb.AF_pbs_overlap.tsv",
                        filter_expr = "CWB_pval <= 0.05 | YHE_pval <= 0.05 | TYD_pval <= 0.05 | WEA_pval <= 0.05")
overlap_and_maybe_write(pbs_gr, Pb_eig1_99_fully, "Pb.AF_fully_pbs_overlap.tsv",
                        filter_expr = "CWB_pval <= 0.05 | YHE_pval <= 0.05 | TYD_pval <= 0.05 | WEA_pval <= 0.05")

# NE with PBS
overlap_and_maybe_write(pbs_gr, NE_eig1_99, "NE.AF_pbs_overlap.tsv",
                        filter_expr = "TYD_pval <= 0.05 | TYE_pval <= 0.05 | WEA_pval <= 0.05")
overlap_and_maybe_write(pbs_gr, NE_eig1_99_fully, "NE.AF_fully_pbs_overlap.tsv",
                        filter_expr = "TYD_pval <= 0.05 | TYE_pval <= 0.05 | WEA_pval <= 0.05")

# Cymru with PBS
overlap_and_maybe_write(pbs_gr, cymru_eig1_99, "cymru.AF_pbs_overlap.tsv",
                        filter_expr = "YHE_pval <= 0.05 | CWB_pval <= 0.05 | AGA_pval <= 0.05")
overlap_and_maybe_write(pbs_gr, cymru_eig1_99_fully, "cymru.AF_fully_pbs_overlap.tsv",
                        filter_expr = "YHE_pval <= 0.05 | CWB_pval <= 0.05 | AGA_pval <= 0.05")

# End of script
