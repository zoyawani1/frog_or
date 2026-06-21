#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tools)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
  stop(
    paste0(
      "\nUsage:\n",
      "  Rscript 06_split_sum_by_hit.R <blast_sum_file> <genome_id> <chunk_base_dir>\n\n",
      "Example:\n",
      "  Rscript 06_split_sum_by_hit.R /path/GCA_xxx_OR_query_1e-20.sum GCA_xxx /path/chunks\n"
    ),
    call. = FALSE
  )
}

blast_sum_file <- args[1]
genome_id      <- args[2]
chunk_base_dir <- args[3]

if (!file.exists(blast_sum_file)) {
  stop(paste("Input blast sum file does not exist:", blast_sum_file), call. = FALSE)
}

chunk_dir <- file.path(chunk_base_dir, genome_id)
dir.create(chunk_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== 06_split_sum_by_hit.R ===\n")
cat("Input sum file:", blast_sum_file, "\n")
cat("Genome ID:", genome_id, "\n")
cat("Chunk output dir:", chunk_dir, "\n")

cat("STEP 1: reading blast sum file\n")
blastsumDD <- read.table(
  blast_sum_file,
  header = TRUE,
  sep = "/",
  quote = "",
  stringsAsFactors = FALSE,
  comment.char = "",
  check.names = FALSE
)

cat("Rows read:", nrow(blastsumDD), "\n")
cat("Cols read:", ncol(blastsumDD), "\n")
cat("Column names:\n")
print(colnames(blastsumDD))
gc()

required_cols <- c(
  "Hit",
  "Hit_length",
  "Percent_identity",
  "Bit score",
  "E-value",
  "Hit_Start",
  "Hit_END",
  "Hit_strand"
)

missing_cols <- setdiff(required_cols, colnames(blastsumDD))
if (length(missing_cols) > 0) {
  stop(
    paste("Missing required columns:", paste(missing_cols, collapse = ", ")),
    call. = FALSE
  )
}

cat("STEP 2: numeric conversion\n")
blastsumDD$Hit_length       <- as.numeric(blastsumDD$Hit_length)
blastsumDD$Percent_identity <- as.numeric(blastsumDD$Percent_identity)
blastsumDD$`Bit score`      <- as.numeric(blastsumDD$`Bit score`)
blastsumDD$`E-value`        <- as.numeric(blastsumDD$`E-value`)
blastsumDD$Hit_Start        <- as.numeric(blastsumDD$Hit_Start)
blastsumDD$Hit_END          <- as.numeric(blastsumDD$Hit_END)
cat("STEP 2 DONE\n")
gc()

cat("STEP 3: filtering\n")
cat("Rows before filtering:", nrow(blastsumDD), "\n")

blastsumDD <- blastsumDD[
  !is.na(blastsumDD$Hit) &
  !is.na(blastsumDD$Hit_length) &
  !is.na(blastsumDD$Percent_identity) &
  !is.na(blastsumDD$`Bit score`) &
  !is.na(blastsumDD$`E-value`) &
  !is.na(blastsumDD$Hit_Start) &
  !is.na(blastsumDD$Hit_END) &
  blastsumDD$Hit_length >= 250 &
  blastsumDD$Percent_identity >= 30 &
  blastsumDD$`Bit score` >= 100 &
  blastsumDD$`E-value` <= 1e-20,
]

cat("Rows after filtering:", nrow(blastsumDD), "\n")
gc()

if (nrow(blastsumDD) == 0) {
  cat("No rows left after filtering. Writing empty marker.\n")
  writeLines(
    genome_id,
    con = file.path(chunk_base_dir, paste0(genome_id, "_NO_ROWS_AFTER_FILTERING.txt"))
  )
  quit(save = "no", status = 0)
}

cat("STEP 4: split by Hit\n")
hits <- unique(blastsumDD$Hit)
cat("Unique hits:", length(hits), "\n")

safe_name <- function(x) {
  gsub("[^A-Za-z0-9._-]", "_", x)
}

manifest_file <- file.path(chunk_dir, paste0(genome_id, "_chunk_manifest.tsv"))
if (file.exists(manifest_file)) file.remove(manifest_file)

manifest_header <- data.frame(
  hit_index    = integer(),
  original_hit = character(),
  safe_hit     = character(),
  chunk_file   = character(),
  n_rows       = integer(),
  stringsAsFactors = FALSE
)

write.table(
  manifest_header,
  file = manifest_file,
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE,
  sep = "\t"
)

for (i in seq_along(hits)) {
  h <- hits[i]
  sub <- blastsumDD[blastsumDD$Hit == h, , drop = FALSE]
  h_safe <- safe_name(h)

  out_file <- file.path(
    chunk_dir,
    paste0(genome_id, "__hit", i, "__", h_safe, ".sum")
  )

  write.table(
    sub,
    file = out_file,
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE,
    sep = "/"
  )

  manifest_row <- data.frame(
    hit_index    = i,
    original_hit = h,
    safe_hit     = h_safe,
    chunk_file   = out_file,
    n_rows       = nrow(sub),
    stringsAsFactors = FALSE
  )

  write.table(
    manifest_row,
    file = manifest_file,
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE,
    sep = "\t",
    append = TRUE
  )
}

cat("Chunk manifest:", manifest_file, "\n")
cat("Total chunk files written:", length(hits), "\n")
cat("=== 06_split_sum_by_hit.R DONE ===\n")
