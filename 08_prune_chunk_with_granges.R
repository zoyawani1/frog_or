#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  require(GenomicRanges, quietly = TRUE)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  stop(
    paste(
      "\n====================",
      "\nUsage:",
      "\nRscript 08_prune_chunk_with_granges.R <chunk_sum_file> <fai_path> <genome_id> <sum_out_dir> <bed_out_dir> [block_size] [boundary_bp]",
      "\n",
      "\narg1 = full path to chunk .sum file",
      "\narg2 = full path to .fai file",
      "\narg3 = genome id",
      "\narg4 = directory for cleaned .sum and info outputs",
      "\narg5 = directory for cleaned .bed outputs",
      "\narg6 = optional row block size (default 50000)",
      "\narg7 = optional boundary bp (default 100)",
      sep = "\n"
    ),
    call. = FALSE
  )
}

chunk_sum_file <- args[1]
fai_path       <- args[2]
genome_id      <- args[3]
sum_out_dir    <- args[4]
bed_out_dir    <- args[5]
block_size     <- if (length(args) >= 6) as.integer(args[6]) else 50000L
boundary_bp    <- if (length(args) >= 7) as.integer(args[7]) else 100L

if (!file.exists(chunk_sum_file)) {
  stop(paste("Input .sum file does not exist:", chunk_sum_file), call. = FALSE)
}

if (!file.exists(fai_path)) {
  stop(paste("FAI file does not exist:", fai_path), call. = FALSE)
}

if (is.na(block_size) || block_size <= 0) {
  stop("block_size must be a positive integer", call. = FALSE)
}

if (is.na(boundary_bp) || boundary_bp < 0) {
  stop("boundary_bp must be a non-negative integer", call. = FALSE)
}

dir.create(sum_out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(bed_out_dir, showWarnings = FALSE, recursive = TRUE)

chunk_label <- sub("\\.sum$", "", basename(chunk_sum_file))

cat("=== BEGIN Rscript 08_prune_chunk_with_granges.R ===\n\n")
cat("Genome ID: ", genome_id, "\n", sep = "")
cat("Input chunk .sum file: ", chunk_sum_file, "\n", sep = "")
cat("Chunk label: ", chunk_label, "\n", sep = "")
cat("FAI path: ", fai_path, "\n", sep = "")
cat("block_size: ", block_size, "\n", sep = "")
cat("boundary_bp: ", boundary_bp, "\n\n", sep = "")

###############################################################################
# STEP 1: read chunk .sum file
###############################################################################
cat("STEP 1: reading chunk .sum file\n")

blastsumDD_raw <- read.table(
  chunk_sum_file,
  header = TRUE,
  sep = "/",
  quote = "",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  comment.char = ""
)

if (nrow(blastsumDD_raw) == 0) {
  cat("No rows in input file: ", chunk_sum_file, "\n", sep = "")
  quit(save = "no")
}

cat("STEP 1 DONE\n")
cat("Rows: ", nrow(blastsumDD_raw), "   Cols: ", ncol(blastsumDD_raw), "\n\n", sep = "")
gc()

###############################################################################
# STEP 2: convert key columns to numeric
###############################################################################
cat("STEP 2: converting key columns to numeric\n")

required_cols <- c(
  "Query",
  "Hit",
  "Hit_length",
  "E-value",
  "Bit score",
  "Percent_identity",
  "Query_Start",
  "Query_End",
  "Hit_Start",
  "Hit_END",
  "Query_strand",
  "Hit_strand"
)

missing_cols <- setdiff(required_cols, colnames(blastsumDD_raw))
if (length(missing_cols) > 0) {
  stop(
    paste("Missing required columns:", paste(missing_cols, collapse = ", ")),
    call. = FALSE
  )
}

blastsumDD_raw$Hit_length       <- suppressWarnings(as.numeric(as.character(blastsumDD_raw$Hit_length)))
blastsumDD_raw$`E-value`        <- suppressWarnings(as.numeric(as.character(blastsumDD_raw$`E-value`)))
blastsumDD_raw$`Bit score`      <- suppressWarnings(as.numeric(as.character(blastsumDD_raw$`Bit score`)))
blastsumDD_raw$Percent_identity <- suppressWarnings(as.numeric(as.character(blastsumDD_raw$Percent_identity)))
blastsumDD_raw$Query_Start      <- suppressWarnings(as.numeric(as.character(blastsumDD_raw$Query_Start)))
blastsumDD_raw$Query_End        <- suppressWarnings(as.numeric(as.character(blastsumDD_raw$Query_End)))
blastsumDD_raw$Hit_Start        <- suppressWarnings(as.numeric(as.character(blastsumDD_raw$Hit_Start)))
blastsumDD_raw$Hit_END          <- suppressWarnings(as.numeric(as.character(blastsumDD_raw$Hit_END)))

cat("STEP 2 DONE\n\n")
gc()

###############################################################################
# STEP 3: length filter
###############################################################################
cat("STEP 3: filtering Hit_length >= 250\n")
cat("Rows before length filter: ", nrow(blastsumDD_raw), "\n", sep = "")

blastsumDD <- blastsumDD_raw[blastsumDD_raw$Hit_length >= 250, , drop = FALSE]

cat("Rows after length filter: ", nrow(blastsumDD), "\n\n", sep = "")
gc()

if (nrow(blastsumDD) == 0) {
  stop("No rows left after Hit_length >= 250 filtering", call. = FALSE)
}

###############################################################################
# STEP 4: NA filtering
###############################################################################
cat("STEP 4: removing rows with missing values in key fields\n")

keep_idx <- !is.na(blastsumDD$Hit_length) &
  !is.na(blastsumDD$`E-value`) &
  !is.na(blastsumDD$`Bit score`) &
  !is.na(blastsumDD$Percent_identity) &
  !is.na(blastsumDD$Hit_Start) &
  !is.na(blastsumDD$Hit_END)

blastsumDD <- blastsumDD[keep_idx, , drop = FALSE]

cat("Rows after NA filter: ", nrow(blastsumDD), "\n\n", sep = "")
gc()

if (nrow(blastsumDD) == 0) {
  cat("No rows left after NA filtering\n")
  quit(save = "no")
}

###############################################################################
# STEP 5: strand normalization and labels
###############################################################################
cat("STEP 5: normalizing strand and creating unique labels\n")

blastsumDD$uniqueHitLabel <- seq_len(nrow(blastsumDD))

blastsumDD$Hit_strand <- trimws(as.character(blastsumDD$Hit_strand))
blastsumDD$Hit_strand <- ifelse(
  blastsumDD$Hit_strand == "1",
  "+",
  ifelse(blastsumDD$Hit_strand == "-1", "-", blastsumDD$Hit_strand)
)

cat("Unique Hit values in this chunk: ", length(unique(blastsumDD$Hit)), "\n\n", sep = "")
gc()

###############################################################################
# STEP 6: prepare row-block chunking inside this single Hit/scaffold
###############################################################################
cat("STEP 6: preparing row-block subchunking\n")

# only used for sorting/block boundaries, not for GRanges itself
blastsumDD$locus_min <- pmin(blastsumDD$Hit_Start, blastsumDD$Hit_END)
blastsumDD$locus_max <- pmax(blastsumDD$Hit_Start, blastsumDD$Hit_END)

ord <- order(blastsumDD$locus_min, blastsumDD$locus_max, na.last = TRUE)
blastsumDD <- blastsumDD[ord, , drop = FALSE]

n_rows <- nrow(blastsumDD)
block_starts <- seq(1L, n_rows, by = block_size)
n_blocks <- length(block_starts)

cat("Rows entering pruning: ", n_rows, "\n", sep = "")
cat("Number of row blocks: ", n_blocks, "\n\n", sep = "")
gc()

###############################################################################
# STEP 7: prune block-by-block using original overlap settings
###############################################################################
cat("STEP 7: pruning block-by-block\n")

cleanByEValue_prune <- function(grobj, overlaps) {
  queryNums <- unique(queryHits(overlaps))
  keep <- GRanges()

  while (length(queryNums) > 0) {
    q <- queryNums[1]
    subs_q <- subjectHits(overlaps[queryHits(overlaps) == q])
    grobj_subs <- grobj[subs_q]
    keepThisHit <- which.min(grobj_subs$evalue)
    keep <- c(keep, grobj_subs[keepThisHit])
    queryNums <- queryNums[!(queryNums %in% subs_q)]
  }

  unique(keep)
}

keep_labels_all <- integer(0)

for (b in seq_along(block_starts)) {
  core_start_idx <- block_starts[b]
  core_end_idx <- min(core_start_idx + block_size - 1L, n_rows)

  core_df <- blastsumDD[core_start_idx:core_end_idx, , drop = FALSE]

  core_min <- min(core_df$locus_min)
  core_max <- max(core_df$locus_max)

  expanded_min <- core_min - boundary_bp
  expanded_max <- core_max + boundary_bp

  sub_idx <- blastsumDD$locus_max >= expanded_min & blastsumDD$locus_min <= expanded_max
  sub_df <- blastsumDD[sub_idx, , drop = FALSE]

  cat("  Block ", b, "/", n_blocks,
      ": core rows ", core_start_idx, "-", core_end_idx,
      " | expanded rows ", nrow(sub_df), "\n", sep = "")

  if (nrow(sub_df) == 0) next

  if (nrow(sub_df) == 1) {
    keep_labels_all <- c(keep_labels_all, sub_df$uniqueHitLabel)
    next
  }

  grobj <- GRanges(
    seqnames = sub_df$Hit,
    ranges = IRanges(
      start = sub_df$Hit_Start,
      end   = sub_df$Hit_END
    ),
    strand = sub_df$Hit_strand,
    bit.score = sub_df$`Bit score`,
    evalue = sub_df$`E-value`,
    query = sub_df$Query,
    hit.length = sub_df$Hit_length,
    uniqueHitLabel = sub_df$uniqueHitLabel
  )

  overlaps <- findOverlaps(
    grobj,
    maxgap = 100L,
    minoverlap = 0L,
    drop.self = FALSE,
    drop.redundant = FALSE
  )

  if (length(queryHits(overlaps)) == 0) {
    kept_sub_labels <- sub_df$uniqueHitLabel
  } else {
    keep_prune_eval <- cleanByEValue_prune(grobj, overlaps)
    kept_sub_labels <- keep_prune_eval$uniqueHitLabel
  }

  # only retain kept labels that belong to the core rows of this block
  kept_core_labels <- intersect(kept_sub_labels, core_df$uniqueHitLabel)
  keep_labels_all <- c(keep_labels_all, kept_core_labels)

  rm(core_df, sub_df, grobj, overlaps, kept_sub_labels, kept_core_labels)
  if (exists("keep_prune_eval")) rm(keep_prune_eval)
  gc()
}

keep_labels_all <- unique(keep_labels_all)

pullOut <- blastsumDD[blastsumDD$uniqueHitLabel %in% keep_labels_all, , drop = FALSE]

cat("\nRows kept after pruning: ", nrow(pullOut), "\n\n", sep = "")

if (nrow(pullOut) == 0) {
  stop("No rows left after pruning", call. = FALSE)
}

gc()

###############################################################################
# STEP 8: prepare BED coordinates
###############################################################################
cat("STEP 8: preparing BED output\n")

pullOut$bed_HitStart0 <- pullOut$Hit_Start - 1
pullOut$bed_HitStart0_minus300 <- pullOut$bed_HitStart0 - 300
pullOut$bed_HitStart0_minus300[pullOut$bed_HitStart0_minus300 < 0] <- 0
pullOut$bed_HitEnd_plus300 <- pullOut$Hit_END + 300

scaffLengths <- read.table(
  fai_path,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE,
  quote = ""
)

if (ncol(scaffLengths) < 2) {
  stop("FAI file must have at least 2 columns", call. = FALSE)
}

colnames(scaffLengths)[1:2] <- c("Scaffold", "length")
scaffLengths <- scaffLengths[, c("Scaffold", "length"), drop = FALSE]
scaffLengths$length <- suppressWarnings(as.numeric(scaffLengths$length))

pullOut_lengths <- merge(
  pullOut,
  scaffLengths,
  by.x = "Hit",
  by.y = "Scaffold",
  all.x = TRUE,
  all.y = FALSE
)

missing_len <- is.na(pullOut_lengths$length)
if (any(missing_len)) {
  stop("Some Hit values in pullOut were not found in the FAI length table", call. = FALSE)
}

over_idx <- pullOut_lengths$bed_HitEnd_plus300 > pullOut_lengths$length
pullOut_lengths$bed_HitEnd_plus300[over_idx] <- pullOut_lengths$Hit_END[over_idx]

pullOut <- pullOut_lengths
rm(pullOut_lengths)
gc()

pullOut$bedstrand <- pullOut$Hit_strand
pullOut$bedname <- paste(
  genome_id,
  pullOut$Hit,
  pullOut$Hit_Start,
  pullOut$Hit_END,
  pullOut$Hit_strand,
  sep = "_"
)

bed <- pullOut[, c("Hit", "bed_HitStart0_minus300", "bed_HitEnd_plus300", "bedname", "E-value", "bedstrand")]
colnames(bed) <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand")

###############################################################################
# STEP 9: write outputs
###############################################################################
cat("STEP 9: writing output files\n")

sum_out  <- file.path(sum_out_dir, paste0(chunk_label, "_cleaned.sum"))
bed_out  <- file.path(bed_out_dir, paste0(chunk_label, "_cleaned.bed"))
info_out <- file.path(sum_out_dir, paste0(chunk_label, "_pullOut.info.txt"))

write.table(
  pullOut,
  file = sum_out,
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE,
  sep = "/"
)
cat("Wrote: ", sum_out, "\n", sep = "")

write.table(
  bed,
  file = bed_out,
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE,
  sep = "\t"
)
cat("Wrote: ", bed_out, "\n", sep = "")

write.table(
  pullOut,
  file = info_out,
  row.names = FALSE,
  col.names = TRUE,
  quote = FALSE,
  sep = "\t"
)
cat("Wrote: ", info_out, "\n", sep = "")

cat("\n=== END Rscript 08_prune_chunk_with_granges.R ===\n")
