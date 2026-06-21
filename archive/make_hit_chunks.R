#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  stop("Usage: Rscript make_hit_chunks.R <fai_file> <sum_file> <out_dir> <n_chunks>", call. = FALSE)
}

fai_file <- args[1]
sum_file <- args[2]
out_dir  <- args[3]
n_chunks <- as.integer(args[4])

if (is.na(n_chunks) || n_chunks < 1) {
  stop("n_chunks must be a positive integer", call. = FALSE)
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== make_hit_chunks.R ===\n")
cat("FAI file: ", fai_file, "\n", sep = "")
cat("SUM file: ", sum_file, "\n", sep = "")
cat("OUT dir : ", out_dir, "\n", sep = "")
cat("Chunks  : ", n_chunks, "\n\n", sep = "")

# Read fai
fai <- read.table(
  fai_file,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = ""
)

if (ncol(fai) < 2) {
  stop("FAI file must have at least 2 columns", call. = FALSE)
}

colnames(fai)[1:2] <- c("Hit", "length")
fai <- fai[, c("Hit", "length")]

# Force length numeric
fai$length <- as.numeric(fai$length)

if (any(is.na(fai$length))) {
  stop("Some scaffold lengths in the FAI file are not numeric", call. = FALSE)
}

# Sort largest to smallest, then greedily balance chunks by total bp
fai <- fai[order(-fai$length), , drop = FALSE]

chunk_bp <- rep(0, n_chunks)
chunk_id <- integer(nrow(fai))

for (i in seq_len(nrow(fai))) {
  j <- which.min(chunk_bp)
  chunk_id[i] <- j
  chunk_bp[j] <- chunk_bp[j] + fai$length[i]
}

fai$chunk <- chunk_id

# Save scaffold -> chunk map
write.table(
  fai,
  file = file.path(out_dir, "chunk_map.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

cat("Chunk bp totals:\n")
for (k in seq_len(n_chunks)) {
  cat("  chunk_", sprintf("%02d", k), ": ", format(chunk_bp[k], big.mark=","), " bp\n", sep = "")
}
cat("\n")

# Read sum file
sum_df <- read.table(
  sum_file,
  header = TRUE,
  sep = "/",
  quote = "",
  stringsAsFactors = FALSE,
  check.names = TRUE,
  comment.char = ""
)

if (!"Hit" %in% colnames(sum_df)) {
  stop("Column 'Hit' not found in .sum file", call. = FALSE)
}

if (nrow(sum_df) == 0) {
  stop("Input .sum file has 0 rows", call. = FALSE)
}

# Preserve row order
sum_df$row_id_tmp__ <- seq_len(nrow(sum_df))

# Merge chunk assignments by Hit
sum_df2 <- merge(
  sum_df,
  fai[, c("Hit", "chunk")],
  by = "Hit",
  all.x = TRUE,
  sort = FALSE
)

# Restore original order
sum_df2 <- sum_df2[order(sum_df2$row_id_tmp__), , drop = FALSE]
sum_df2$row_id_tmp__ <- NULL

# Check for Hit names missing from fai
if (any(is.na(sum_df2$chunk))) {
  missing_hits <- unique(sum_df2$Hit[is.na(sum_df2$chunk)])
  writeLines(missing_hits, file.path(out_dir, "missing_hits_from_fai.txt"))
  stop("Some Hit values in the .sum file were not found in the FAI. See missing_hits_from_fai.txt", call. = FALSE)
}

# Write chunk sum files
for (k in seq_len(n_chunks)) {
  df_k <- sum_df2[sum_df2$chunk == k, setdiff(colnames(sum_df2), "chunk"), drop = FALSE]
  out_file <- file.path(out_dir, sprintf("chunk_%02d.sum", k))

  write.table(
    df_k,
    file = out_file,
    sep = "/",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE
  )

  cat("Wrote ", out_file, " with ", nrow(df_k), " rows\n", sep = "")
}

cat("\nDone.\n")
