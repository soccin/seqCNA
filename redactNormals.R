require(tidyverse)

#
# SeqCNA labels each sample as "<tumor>__<matched normal>". Strip the
# "__<matched normal>" half from the project output files so the normals
# are not disclosed. Every file is rewritten in place.
#
# write_xlsx() is the openxlsx wrapper defined in ~/.Rprofile
#

normal_id_pattern <- "__.*"

#' Locate the one project output file matching a pattern
#'
#' @param pattern regexp matched against paths below the working directory
#' @return path of the matching file
find_output_file <- function(pattern) {
  files <- fs::dir_ls(".", regexp = pattern, recurse = TRUE)

  # Rewriting in place, so anything other than a single hit is unsafe
  if (length(files) != 1) {
    rlang::abort(
      str_glue("Expected 1 file matching '{pattern}', found {length(files)}")
    )
  }

  files
}

#
# These files name the sample in an ID column
#

gene_table_file <- find_output_file("___GeneTable.xlsx")

gene_table_file |>
  readxl::read_xlsx() |>
  mutate(ID = str_remove(ID, normal_id_pattern)) |>
  write_xlsx(gene_table_file)

igv_seg_file <- find_output_file("___IGV.seg")

igv_seg_file |>
  read_tsv(show_col_types = FALSE) |>
  mutate(ID = str_remove(ID, normal_id_pattern)) |>
  write_tsv(igv_seg_file)

#
# These files name the samples in the column headers
#

gene_matrix_file <- find_output_file("GeneMatrix.csv")

gene_matrix_file |>
  read_csv(show_col_types = FALSE) |>
  rename_with(~ str_remove(., normal_id_pattern)) |>
  write_csv(gene_matrix_file)

segment_matrix_file <- find_output_file("SegmentMatrix.csv")

segment_matrix_file |>
  read_csv(show_col_types = FALSE) |>
  rename_with(~ str_remove(., normal_id_pattern)) |>
  write_csv(segment_matrix_file)
