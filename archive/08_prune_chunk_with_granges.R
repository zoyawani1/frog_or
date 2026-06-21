#!/usr/bin/env Rscript

###############################################################################
# Step 2: prune one chunked .sum file (one Hit/scaffold per file)
#
# Logic is based on the original script, with only these changes:
# - input is a single chunk .sum file instead of one full-genome .sum
# - full .fai path is passed directly
# - outputs go to step-2 output directories
# - key numeric columns are explicitly converted earlier for safety
# - safer BED overhang reassignment is used
#
# Kept from original logic:
# - use Hit_Start and Hit_END directly inside IRanges
# - use original overlap settings:
#     maxgap=100L, minoverlap=0L, drop.self=FALSE, drop.redundant=FALSE
# - prune by minimum E-value only
# - use original strand normalization idea
# - no scaffold loop
###############################################################################

suppressPackageStartupMessages({
  require(GenomicRanges, quietly = TRUE)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  stop(
    paste(
      "\n====================",
      "\nUsage:",
      "\nRscript 08_prune_chunk_with_granges.R <chunk_sum_file> <fai_path> <genome_id> <sum_out_dir> <bed_out_dir>",
      "\n",
      "\narg1 = full path to chunk .sum file",
      "\narg2 = full path to .fai file",
      "\narg3 = genome id",
      "\narg4 = directory for cleaned .sum and info outputs",
      "\narg5 = directory for cleaned .bed outputs",
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

if (!file.exists(chunk_sum_file)) {
  stop(paste("Input .sum file does not exist:", chunk_sum_file), call. = FALSE)
}

if (!file.exists(fai_path)) {
  stop(paste("FAI file does not exist:", fai_path), call. = FALSE)
}

dir.create(sum_out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(bed_out_dir, showWarnings = FALSE, recursive = TRUE)

chunk_label <- sub("\\.sum$", "", basename(chunk_sum_file))

cat("=== BEGIN Rscript 08_prune_chunk_with_granges.R ===\n\n")
cat("Genome ID: ", genome_id, "\n", sep = "")
cat("Input chunk .sum file: ", chunk_sum_file, "\n", sep = "")
cat("Chunk label: ", chunk_label, "\n", sep = "")
cat("FAI path: ", fai_path, "\n\n", sep = "")

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
# STEP 3: length filter (repeat defensively, even though step 1 already did it)
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
# STEP 4: NA filtering (kept after length filter to match original order)
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
# STEP 5: strand normalization and row labels
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
# STEP 6: build GRanges using original coordinate logic
# (use Hit_Start and Hit_END directly, no pmin/pmax normalization)
###############################################################################
cat("STEP 6: building GRanges\n")

grobj <- GRanges(
  seqnames = blastsumDD$Hit,
  ranges = IRanges(
    start = blastsumDD$Hit_Start,
    end   = blastsumDD$Hit_END
  ),
  strand = blastsumDD$Hit_strand,
  bit.score = blastsumDD$`Bit score`,
  evalue = blastsumDD$`E-value`,
  query = blastsumDD$Query,
  hit.length = blastsumDD$Hit_length,
  uniqueHitLabel = blastsumDD$uniqueHitLabel
)

cat("STEP 6 DONE\n\n")
gc()

###############################################################################
# STEP 7: find overlaps using original settings
###############################################################################
cat("STEP 7: finding overlaps\n")

overlaps <- findOverlaps(
  grobj,
  maxgap = 100L,
  minoverlap = 0L,
  drop.self = FALSE,
  drop.redundant = FALSE
)

cat("Number of overlap pairs: ", length(queryHits(overlaps)), "\n\n", sep = "")

###############################################################################
# STEP 8: prune by minimum E-value only
###############################################################################
cat("STEP 8: pruning by minimum E-value\n")

cleanByEValue_prune <- function(grobj, overlaps) {
  queryNums <- unique(queryHits(overlaps))
  keep <- GRanges()

  total_queries <- length(queryNums)
  if (total_queries > 0) {
    pb <- txtProgressBar(min = 0, max = total_queries, style = 3)
  }

  while (length(queryNums) > 0) {
    q <- queryNums[1]

    # original logic keeps hits-to-self inside overlaps
    subs_q <- subjectHits(overlaps[queryHits(overlaps) == q])
    grobj_subs <- grobj[subs_q]

    keepThisHit <- which.min(grobj_subs$evalue)
    keep <- c(keep, grobj_subs[keepThisHit])

    queryNums <- queryNums[!(queryNums %in% subs_q)]

    if (total_queries > 0) {
      setTxtProgressBar(pb, total_queries - length(queryNums))
    }
  }

  unique(keep)
}

if (length(queryHits(overlaps)) == 0) {
  keep_labels <- blastsumDD$uniqueHitLabel
} else {
  keep_prune_eval <- cleanByEValue_prune(grobj, overlaps)
  keep_labels <- keep_prune_eval$uniqueHitLabel
}

pullOut <- blastsumDD[blastsumDD$uniqueHitLabel %in% keep_labels, , drop = FALSE]

cat("\nRows kept after pruning: ", nrow(pullOut), "\n\n", sep = "")

if (nrow(pullOut) == 0) {
  stop("No rows left after pruning", call. = FALSE)
}

gc()

###############################################################################
# STEP 9: prepare BED coordinates
###############################################################################
cat("STEP 9: preparing BED output\n")

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

# safer equivalent of the original overhang fix
over_idx <- pullOut_lengths$bed_HitEnd_plus300 > pullOut_lengths$length
pullOut_lengths$bed_HitEnd_plus300[over_idx] <- pullOut_lengths$Hit_END[over_idx]

pullOut <- pullOut_lengths
rm(pullOut_lengths)
gc()

pullOut$bedstrand <- pullOut$Hit_strand

# best for later merge-per-genome step:
# keep row names informative but not dependent on chunk filename
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
# STEP 10: write outputs
###############################################################################
cat("STEP 10: writing output files\n")

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
