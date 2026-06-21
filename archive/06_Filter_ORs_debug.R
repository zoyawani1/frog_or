#!/usr/bin/env Rscript

###############################################################################
# Filter OR BLAST hits for one genome or one chunk of a genome.
# Input .sum is slash-delimited and must contain a column named "Hit".
# "Hit" is treated as the scaffold/chromosome/contig identifier.
###############################################################################

suppressPackageStartupMessages({
  require(GenomicRanges, quietly = TRUE)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop(paste(
    "\n====================",
    "\nUsage:",
    "\nRscript 06_Filter_ORs_debug.R <SPP_PATH> <fai_name> [blast_sum_file] [chunk_label]",
    "\n",
    "\nargs1 - absolute/path/to/Species",
    "\nargs2 - fna.fai filename located inside SPP_PATH",
    "\nargs3 - optional full path to input .sum file",
    "\nargs4 - optional chunk label",
    sep = "\n"
  ), call. = FALSE)
}

SPP_PATH <- args[1]
fai_name <- args[2]
genome_id <- basename(SPP_PATH)

BLAST_DIR <- "/project/stuckert/users/wani/f2/anuran_genomes/blast_results_1e-20"

if (length(args) >= 3 && nzchar(args[3])) {
  blast_sum_file <- args[3]
} else {
  blast_sum_file <- file.path(BLAST_DIR, paste0(genome_id, "_OR_query_1e-20.sum"))
}

if (length(args) >= 4 && nzchar(args[4])) {
  chunk_label <- args[4]
} else {
  chunk_label <- "full"
}

BED_BASE_DIR <- "/project/stuckert/users/wani/f2/anuran_genomes/bed_files"
dir.create(BED_BASE_DIR, showWarnings = FALSE, recursive = TRUE)

setwd(SPP_PATH)

cat("=== BEGIN Rscript Filter_ORs.R ===\n\n")
cat("Genome ID: ", genome_id, "\n", sep = "")
cat("Input sum file: ", blast_sum_file, "\n", sep = "")
cat("Chunk label: ", chunk_label, "\n\n", sep = "")

if (!file.exists(blast_sum_file)) {
  stop(paste("Input .sum file does not exist:", blast_sum_file), call. = FALSE)
}

fai_path <- file.path(SPP_PATH, fai_name)
if (!file.exists(fai_path)) {
  stop(paste("FAI file does not exist:", fai_path), call. = FALSE)
}

###############################################################################
# STEP 1: Read .sum file
###############################################################################
cat("STEP 1: reading .sum file\n")

blastsumDD_raw <- read.table(
  blast_sum_file,
  header = TRUE,
  sep = "/",
  quote = "",
  stringsAsFactors = FALSE,
  check.names = TRUE,
  comment.char = ""
)

if (nrow(blastsumDD_raw) == 0) {
  cat("No rows in input file: ", blast_sum_file, "\n", sep = "")
  quit(save = "no")
}

cat("STEP 1 DONE\n")
cat("Rows: ", nrow(blastsumDD_raw), "  Cols: ", ncol(blastsumDD_raw), "\n\n", sep = "")
gc()

###############################################################################
# STEP 2: convert key columns to numeric
###############################################################################
cat("STEP 2: converting numeric columns\n")

required_cols <- c(
  "Hit",
  "Hit_length",
  "E.value",
  "Bit.score",
  "Percent_identity",
  "Query_Start",
  "Query_End",
  "Hit_Start",
  "Hit_END",
  "Hit_strand"
)

missing_cols <- setdiff(required_cols, colnames(blastsumDD_raw))
if (length(missing_cols) > 0) {
  stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")), call. = FALSE)
}

blastsumDD_raw$Hit_length <- as.numeric(as.character(blastsumDD_raw$Hit_length))
blastsumDD_raw$Percent_identity <- as.numeric(as.character(blastsumDD_raw$Percent_identity))
blastsumDD_raw$E.value <- as.numeric(as.character(blastsumDD_raw$E.value))
blastsumDD_raw$Bit.score <- as.numeric(as.character(blastsumDD_raw$Bit.score))
blastsumDD_raw$Query_Start <- as.numeric(as.character(blastsumDD_raw$Query_Start))
blastsumDD_raw$Query_End <- as.numeric(as.character(blastsumDD_raw$Query_End))
blastsumDD_raw$Hit_Start <- as.numeric(as.character(blastsumDD_raw$Hit_Start))
blastsumDD_raw$Hit_END <- as.numeric(as.character(blastsumDD_raw$Hit_END))
blastsumDD_raw$Hit_strand <- as.numeric(as.character(blastsumDD_raw$Hit_strand))

cat("STEP 2 DONE\n\n")
gc()



###############################################################################
# STEP 3: remove NA rows in key fields
###############################################################################
cat("STEP 3: removing rows with missing values in key fields\n")

keep_idx <- !is.na(blastsumDD$Hit_length) &
            !is.na(blastsumDD$Percent_identity) &
            !is.na(blastsumDD$E.value) &
            !is.na(blastsumDD$Bit.score) &
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
# STEP 4: normalize coordinates and strand
###############################################################################
cat("STEP 4: normalizing coordinates and strand\n")

blastsumDD$uniqueHitLabel <- seq_len(nrow(blastsumDD))

blastsumDD$Hit_strand <- ifelse(
  blastsumDD$Hit_strand == 1, "+",
  ifelse(blastsumDD$Hit_strand == -1, "-", as.character(blastsumDD$Hit_strand))
)

blastsumDD$locus_start <- pmin(blastsumDD$Hit_Start, blastsumDD$Hit_END)
blastsumDD$locus_end   <- pmax(blastsumDD$Hit_Start, blastsumDD$Hit_END)

cat("STEP 5 DONE\n\n")
gc()

###############################################################################
# STEP 5: quality filtering
###############################################################################
cat("STEP 5: filtering by E-value, hit length, and percent identity\n")
cat("Rows before quality filtering: ", nrow(blastsumDD), "\n", sep = "")

blastsumDD <- blastsumDD[
  !is.na(blastsumDD$Hit_length) &
  !is.na(blastsumDD$Percent_identity) &
  !is.na(blastsumDD$E.value) &
  !is.na(blastsumDD$Bit.score) &
  blastsumDD$Hit_length >= 250 &
  blastsumDD$Percent_identity >= 30 &
  blastsumDD$E.value <= 1e-20,
  ,
  drop = FALSE
]

cat("Rows after quality filtering: ", nrow(blastsumDD), "\n\n", sep = "")

if (nrow(blastsumDD) == 0) {
  stop("No rows left after Step 6.5 quality filtering", call. = FALSE)
}

gc()


###############################################################################
# STEP 6: scaffold-by-scaffold pruning with GRanges
###############################################################################
cat("STEP 6: scaffold-by-scaffold pruning\n")
cat("Number of rows going into pruning: ", nrow(blastsumDD), "\n", sep = "")

cleanByEValue_prune <- function(grObj, overlaps) {
  queryNums <- unique(queryHits(overlaps))
  keep <- GRanges()

  while (length(queryNums) > 0) {
    q <- queryNums[1]
    subs_q <- subjectHits(overlaps[queryHits(overlaps) == q])
    grobj_subs <- grObj[subs_q]
    keepThisHit <- which.min(grobj_subs$evalue)
    keep <- c(keep, grobj_subs[keepThisHit])
    queryNums <- queryNums[!(queryNums %in% subs_q)]
  }

  unique(keep)
}

scaffolds <- unique(blastsumDD$Hit)
cat("Number of unique scaffolds: ", length(scaffolds), "\n", sep = "")

keep_labels <- character(0)

for (scaf in scaffolds) {
  cat("Processing scaffold: ", scaf, "\n", sep = "")

  scaf_df <- blastsumDD[blastsumDD$Hit == scaf, , drop = FALSE]
  cat("Rows on scaffold: ", nrow(scaf_df), "\n", sep = "")

  if (nrow(scaf_df) == 0) next

  if (nrow(scaf_df) == 1) {
    keep_labels <- c(keep_labels, scaf_df$uniqueHitLabel)
    next
  }

  grobj_scaf <- GRanges(
    seqnames = scaf_df$Hit,
    ranges = IRanges(start = scaf_df$Hit_Start, end = scaf_df$Hit_END),
    strand = scaf_df$Hit_strand,
    bit.score = scaf_df$Bit.score,
    evalue = scaf_df$E.value,
    query = scaf_df$Query,
    hit.length = scaf_df$Hit_length,
    uniqueHitLabel = scaf_df$uniqueHitLabel
  )

  overlaps_scaf <- findOverlaps(
    grobj_scaf,
    maxgap = 100L,
    minoverlap = 0L,
    drop.self = TRUE,
    drop.redundant = TRUE
  )

  cat("Overlap pairs on scaffold: ", length(queryHits(overlaps_scaf)), "\n", sep = "")

  if (length(queryHits(overlaps_scaf)) == 0) {
    keep_labels <- c(keep_labels, scaf_df$uniqueHitLabel)
    rm(scaf_df, grobj_scaf, overlaps_scaf)
    gc()
    next
  }

  keep_prune_eval <- cleanByEValue_prune(grobj_scaf, overlaps_scaf)
  keep_labels <- c(keep_labels, keep_prune_eval$uniqueHitLabel)

  rm(scaf_df, grobj_scaf, overlaps_scaf, keep_prune_eval)
  gc()
}

pullOut <- blastsumDD[blastsumDD$uniqueHitLabel %in% unique(keep_labels), , drop = FALSE]

cat("STEP 7 DONE\n")
cat("Number kept after pruning: ", nrow(pullOut), "\n\n", sep = "")

if (nrow(pullOut) == 0) {
  stop("No rows left after scaffold pruning", call. = FALSE)
}

###############################################################################
# STEP 8: prepare BED coordinates
###############################################################################
cat("STEP 8: preparing BED output\n")

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
scaffLengths <- scaffLengths[, c("Scaffold", "length")]
scaffLengths$length <- as.numeric(scaffLengths$length)

pullOut$bed_HitStart0 <- pullOut$Hit_Start - 1
pullOut$bed_HitStart0_minus300 <- pullOut$bed_HitStart0 - 300
pullOut$bed_HitStart0_minus300[pullOut$bed_HitStart0_minus300 < 0] <- 0
pullOut$bed_HitEnd_plus300 <- pullOut$Hit_END + 300

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
  pullOut$Hit,
  pullOut$Hit_Start,
  pullOut$Hit_END,
  pullOut$Hit_strand,
  sep = "_"
)

bed <- pullOut[, c("Hit", "bed_HitStart0_minus300", "bed_HitEnd_plus300", "bedname", "E.value", "bedstrand")]
colnames(bed) <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand")

###############################################################################
# STEP 9: write outputs
###############################################################################
cat("STEP 9: writing output files\n")

sum_out <- file.path(BED_BASE_DIR, paste0(genome_id, "_", chunk_label, "_cleaned.sum"))
bed_out <- file.path(BED_BASE_DIR, paste0(genome_id, "_", chunk_label, "_cleaned.bed"))
info_out <- file.path(BED_BASE_DIR, paste0(genome_id, "_", chunk_label, "_pullout.info.txt"))

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

cat("\n=== END Rscript Filter_ORs.R ===\n")
