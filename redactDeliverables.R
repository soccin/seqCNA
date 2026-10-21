suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

#
# SeqCNA labels each sample by its pair id "<tumor>__<normal>". Replace the
# pair id by the tumor id in the deliverables of PROJDIR, so they do not
# show which normal was used. The map is PAIRING (pairing.tsv, written by
# selectBestMatch: TUMOR, NORMAL, PAIR, one row per tumor), not a regex: a
# sample can be both a tumor and a normal. Rewrites in place
#
#   <ProjNo>___IGV.seg          ID column
#   <ProjNo>___GeneTable.xlsx   ID column
#   GeneMatrix.csv              sample column headers
#   SegmentMatrix.csv           sample column headers
#
# The xlsx and the CSVs exist only when getGeneCalls found a significant
# gene (otherwise it writes README.txt). The PNGs keep the pair id.
#
# Run by postProcess.sh with: bin/rscript.sh redactDeliverables.R PAIRING PROJDIR
#

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  cat("\n   usage: seqCNA/bin/rscript.sh seqCNA/redactDeliverables.R PAIRING PROJDIR\n\n")
  quit(status = 1)
}
pairing_file <- args[1]
proj_dir <- args[2]
proj_no <- fs::path_file(proj_dir)

pairing <- read_tsv(pairing_file, col_types = cols(.default = "c"))

#
# One best normal per tumor, so the tumor ids are unique; otherwise two
# samples would be merged under one name
#
dup_tumors <- pairing$TUMOR[duplicated(pairing$TUMOR)]
if (length(dup_tumors) > 0 || anyDuplicated(pairing$PAIR) > 0) {
  rlang::abort(str_glue(
    "{pairing_file}: tumor or pair listed twice: {str_c(dup_tumors, collapse = ' ')}"
  ))
}

pair_to_tumor <- rlang::set_names(pairing$TUMOR, pairing$PAIR)

#' Map pair ids to tumor ids; every id must be in the pairing file
#'
#' @param ids pair ids
#' @param where file the ids come from, for the error message
#' @return tumor ids, in the same order
redact_ids <- function(ids, where) {
  unknown <- setdiff(ids, names(pair_to_tumor))
  if (length(unknown) > 0) {
    rlang::abort(str_glue(
      "{where}: not in {pairing_file}: {str_c(unknown, collapse = ' ')}"
    ))
  }
  unname(pair_to_tumor[ids])
}

#
# Text files are read and written as character so that the numbers are
# written back exactly as they were
#
read_as_text <- function(file, delim) {
  read_delim(file, delim = delim, col_types = cols(.default = "c"),
             na = character(), trim_ws = FALSE)
}

#
# These files name the sample in an ID column
#

igv_seg_file <- fs::path(proj_dir, str_glue("{proj_no}___IGV.seg"))

igv_seg_file |>
  read_as_text("\t") |>
  mutate(ID = redact_ids(ID, igv_seg_file)) |>
  write_tsv(igv_seg_file, na = "")

gene_table_file <- fs::path(proj_dir, str_glue("{proj_no}___GeneTable.xlsx"))

if (fs::file_exists(gene_table_file)) {
  gene_table_file |>
    readxl::read_xlsx() |>
    mutate(ID = redact_ids(ID, gene_table_file)) |>
    openxlsx::write.xlsx(gene_table_file)
}

#
# These files name the samples in the column headers; getGeneCalls names
# the sample columns by pair id, which starts with "s_"
#

redact_columns <- function(file) {
  if (!fs::file_exists(file)) {
    return(invisible(NULL))
  }
  file |>
    read_as_text(",") |>
    rename_with(~ redact_ids(.x, file), matches("^s_")) |>
    write_csv(file, na = "")
}

redact_columns(fs::path(proj_dir, "GeneMatrix.csv"))
redact_columns(fs::path(proj_dir, "SegmentMatrix.csv"))
