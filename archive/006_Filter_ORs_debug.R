#!/usr/bin/env Rscript

#################################################################################
# WRITE FUNCTION COMMENTS HERE 
#################################################################################

require(GenomicRanges, quietly = TRUE)
# try(library(GenomicRanges), silent=TRUE)

args <- commandArgs(trailingOnly = TRUE)

# Check for if arguments are specified. 
if (length(args) < 2) {
    stop("\n=================\nTwo must be supplied. \n
         args1 - absolute/path/to/Species. \n
         args2 - fna.fai filename. (file should be in Species directory above. \nGenerated via 1_ab_tblastN_ORs. )", call.=FALSE)
}  else if(!file.exists(file.path( args[1], args[2]))){
    stop("\n=================\nThe provided argument file path of ", 
         file.path( args[1], args[2]),
         " does not exist.\n Please provide a fai filename in ", args[1],call.=FALSE)
}

#NEW 3/25
SPP_PATH <- args[1]
fai_name <- args[2]
spp <- basename(SPP_PATH)

BED_BASE_DIR <- "/project/stuckert/users/wani/f2/anuran_genomes/bed_files"
dir.create(BED_BASE_DIR, showWarnings = FALSE, recursive = TRUE)

setwd(SPP_PATH)

cat("=== BEGIN Rscript Filter_ORs.R ===\n\n")




#original
# Set arguments to variables 
#SPP_PATH=args[1]
#setwd(SPP_PATH)
#spp=tail(unlist(strsplit(SPP_PATH, '/')), 1) # Grabs the last directory as spp, should be Specie

#cat("=== BEGIN Rscript Filter_ORs.R ===\n\n")

########### read back in the deduped (DD) results that you just wrote out (once) #########
# Read in result generated from 1_ab_tblastN_OR.pl
#cat(paste0("=== Read in result generated from 1_ab_tblastN_OR.pl === \nfilepath: ",file.path( SPP_PATH, "OR_query_mini_ORN_gharial.sum")))
#blastsumDD_raw <- read.table(file.path( SPP_PATH, "OR_query_mini_ORN_gharial.sum"),header=T,sep="/",quote="")
########### read back in the deduped (DD) results ###########
##New section added becuase original script hard codes the sum file 
#sum_files <- list.files(
#  path = SPP_PATH,
#  pattern = "_OR_query\\.sum$",
#  full.names = TRUE
#)
#
#if (length(sum_files) == 0) {
#  stop("No *_OR_query.sum file found in ", SPP_PATH, call.=FALSE)
#}

#if (length(sum_files) > 1) {
#  stop(
#    "Multiple *_OR_query.sum files found:\n",
#    paste(basename(sum_files), collapse = "\n"),
#    "\nPlease keep only one per genome directory.",
#    call.=FALSE
#  )
#}

#blast_sum_file <- sum_files[1]

######################NEW SECTION BECAUSE ORIGINAL SCRIPT CAN ONLY HANDLE ONE SUM FILE, INSTEAD OF LOOPING THROUGH MULTIPLE##


BLAST_DIR <- "/project/stuckert/users/wani/f2/anuran_genomes/blast_results_1e-20"

genome_id <- basename(SPP_PATH)

blast_sum_file <- file.path(BLAST_DIR, paste0(genome_id, "_OR_query_1e-20.sum"))

if (!file.exists(blast_sum_file)) {

    cat("WARNING: BLAST sum file not found for genome:", genome_id, "\n")

    write(
        genome_id,
        file = "/project/stuckert/users/wani/f2/anuran_genomes/bed_files/logs/missing_sum_files.txt",
        append = TRUE
    )

    quit(save="no")
}
################# NEW ADDITION TO FIGURE OUT WHICH STEP CAUSES OOM ERR ##################################
cat("=== Using BLAST sum file ===\n", blast_sum_file, "\n")
cat("STEP 1: starting read.table\n")

blastsumDD_raw <- read.table(
    blast_sum_file,
    header = TRUE,
    sep = "/",
    quote = ""
)

cat("STEP 1 DONE\n")
cat("Rows:", nrow(blastsumDD_raw), "Cols:", ncol(blastsumDD_raw), "\n")
gc()

cat("STEP 2: converting Hit_length\n")
blastsumDD_raw$Hit_length <- as.numeric(as.character(blastsumDD_raw$Hit_length))
cat("STEP 2 DONE\n")
gc()

cat("STEP 3: filtering by length >= 250\n")
blastsumDD <- blastsumDD_raw[blastsumDD_raw$Hit_length >= 250,]
cat("STEP 3 DONE\n")
cat("Rows after filter:", nrow(blastsumDD), "\n")
gc()

cat("STEP 4: removing NA\n")
blastsumDD <- na.omit(blastsumDD)
cat("STEP 4 DONE\n")
cat("Rows after na.omit:", nrow(blastsumDD), "\n")
gc()
cat("STEP 5: converting coordinates\n")
blastsumDD$Hit_Start <- as.numeric(as.character(blastsumDD$Hit_Start))
blastsumDD$Hit_END <- as.numeric(as.character(blastsumDD$Hit_END))
cat("STEP 5 DONE\n")

rm(blastsumDD_raw)
gc()
##########################################################################################################


#ORIGINAL
#cat("=== Using BLAST sum file ===\n", blast_sum_file, "\n")

###blastsumDD_raw <- read.table(
#  blast_sum_file,
#  header = TRUE,
#  sep = "/",
#  quote = ""
#)

# do some filtering ahead of time to speed up processing:
# get rids of hits < 250bp in length
# tail(blastsumDD_raw)
# dim(blastsumDD_raw)

#blastsumDD_raw$Hit_length <- as.numeric(as.character(blastsumDD_raw$Hit_length))
#blastsumDD <- blastsumDD_raw[blastsumDD_raw$Hit_length >= 250,]
#blastsumDD <- na.omit(blastsumDD)
#blastsumDD$Hit_Start <- as.numeric(as.character(blastsumDD$Hit_Start))
#blastsumDD$Hit_END <- as.numeric(as.character(blastsumDD$Hit_END))

# rownames(blastsumDD) <- NULL
blastsumDD$uniqueHitLabel <- seq(1:dim(blastsumDD)[1])

# Converts strands to '+' and '-'
blastsumDD$Hit_strand <- ifelse(blastsumDD$Hit_strand == '1', 
                                '+', 
                                ifelse(blastsumDD$Hit_strand == '-1', 
                                        '-', blastsumDD$Hit_strand))
cat("STEP 6  DONE\n")
#################################################################################
# 1. need to remove overlapping hits (within 100bp of each other) by choosing 
# the one with the best e-value (lowest) or bit score (highest) tail(blastsumDD) 
# so what matters is Hit (scaffold ID), Hit_Start, Hit_End, Evalue and Bit Score of each hit.

# source('https://bioconductor.org/biocLite.R')
# biocLite('GenomicRanges')

# seqnames are chromosomes/scaffolds
# ranges use IRanges to give starts then ends 
# strand you get from results
# score is bit score (or could be evalue)
#################################################################################

#write.csv(as.data.frame(blastsumDD), file="blastsumDD_test.csv")
###############ORIGINAL CODE #######################################

#cat("STEP 7: creating GRanges\n")
#cat("Number of rows going into pruning:", nrow(blastsumDD), "\n")
#gc()
#
#grobj <- GRanges(seqnames=blastsumDD$Hit,
#                 ranges=IRanges(start=blastsumDD$Hit_Start,
#                                end=blastsumDD$Hit_END),
#                 strand=blastsumDD$Hit_strand,
#                 bit.score=blastsumDD$Bit.score,
#                 evalue=blastsumDD$E.value,
#                query=blastsumDD$Query,
#                 hit.length=blastsumDD$Hit_length,
#                 uniqueHitLabel=blastsumDD$uniqueHitLabel)
#
#cat("STEP 7  DONE\n")
#####################################################NEW PRUNING FILTER 3/31/2026###########################
cat("STEP 6.5: filtering by E-value, hit length, and percent identity\n")

blastsumDD$Hit_length <- as.numeric(as.character(blastsumDD$Hit_length))
blastsumDD$Percent_identity <- as.numeric(as.character(blastsumDD$Percent_identity))
blastsumDD$E.value <- as.numeric(as.character(blastsumDD$E.value))
blastsumDD$Bit.score <- as.numeric(as.character(blastsumDD$Bit.score))

cat("Rows before quality filtering:", nrow(blastsumDD), "\n")

blastsumDD <- blastsumDD[
  !is.na(blastsumDD$Hit_length) &
  !is.na(blastsumDD$Percent_identity) &
  !is.na(blastsumDD$E.value) &
  !is.na(blastsumDD$Bit.score) &
  blastsumDD$Hit_length >= 250 &
  blastsumDD$Percent_identity >= 30 &
  blastsumDD$E.value <= 1e-20,
]

cat("Rows after quality filtering:", nrow(blastsumDD), "\n")

if (nrow(blastsumDD) == 0) {
    stop("No rows left after Step 6.5 quality filtering", call. = FALSE)
}

gc()

#################################NEW PRUNING 4/8#################################
cat("STEP 6.6: defining loci and keeping top 10 hits per locus\n")

# Make sure coordinate columns are numeric
blastsumDD$Hit_Start <- as.numeric(as.character(blastsumDD$Hit_Start))
blastsumDD$Hit_END <- as.numeric(as.character(blastsumDD$Hit_END))

# Define locus coordinates robustly in case strand flips start/end
blastsumDD$locus_start <- pmin(blastsumDD$Hit_Start, blastsumDD$Hit_END)
blastsumDD$locus_end   <- pmax(blastsumDD$Hit_Start, blastsumDD$Hit_END)

# Sort by scaffold and genomic position
blastsumDD <- blastsumDD[order(blastsumDD$Hit, blastsumDD$locus_start, blastsumDD$locus_end), ]

# Distance allowed between neighboring hits in the same locus
LOCUS_GAP <- 10000

assign_loci <- function(df, gap = 10000) {
    if (nrow(df) == 0) return(df)

    locus_id <- integer(nrow(df))
    current_locus <- 1
    current_end <- df$locus_end[1]
    locus_id[1] <- current_locus

    if (nrow(df) > 1) {
        for (i in 2:nrow(df)) {
            if (df$locus_start[i] > (current_end + gap)) {
                current_locus <- current_locus + 1
                current_end <- df$locus_end[i]
            } else {
                current_end <- max(current_end, df$locus_end[i])
            }
            locus_id[i] <- current_locus
        }
    }

    df$locus_id <- paste(df$Hit[1], locus_id, sep = "_")
    df
}

blastsumDD_loci <- do.call(
    rbind,
    lapply(split(blastsumDD, blastsumDD$Hit), assign_loci, gap = LOCUS_GAP)
)

rownames(blastsumDD_loci) <- NULL

cat("Rows before top-10-per-locus filter:", nrow(blastsumDD_loci), "\n")
cat("Number of loci before top-10-per-locus filter:", length(unique(blastsumDD_loci$locus_id)), "\n")

blastsumDD_top10 <- do.call(
    rbind,
    lapply(split(blastsumDD_loci, blastsumDD_loci$locus_id), function(df) {
        df <- df[order(
            -df$Bit.score,
            df$E.value,
            -df$Percent_identity,
            -df$Hit_length
        ), ]
        head(df, 10)
    })
)

rownames(blastsumDD_top10) <- NULL

blastsumDD <- blastsumDD_top10

# Reset uniqueHitLabel AFTER filtering
blastsumDD$uniqueHitLabel <- seq_len(nrow(blastsumDD))

cat("Rows after top-10-per-locus filter:", nrow(blastsumDD), "\n")
cat("Number of loci after top-10-per-locus filter:", length(unique(blastsumDD$locus_id)), "\n")
cat("Number of unique scaffolds after top-10-per-locus filter:", length(unique(blastsumDD$Hit)), "\n")


if (nrow(blastsumDD) == 0) {
    stop("No rows left after Step 6.6 top-10-per-locus filtering", call. = FALSE)
}

gc()
cat("STEP 6.6 DONE\n")

#################################################END OF NEW PRUNING 4/8###############################################




###################NEW PRUNING 4/2#######################
#cat("STEP 6.6: keeping top 10 hits per scaffold by Bit score\n")
#
#blastsumDD <- blastsumDD[order(
#    blastsumDD$Hit,
#    -blastsumDD$Bit.score
#), ]
#
#blastsumDD <- do.call(rbind, lapply(split(blastsumDD, blastsumDD$Hit), function(df) {
#  head(df, 10)
#}))
#
#rownames(blastsumDD) <- NULL
#blastsumDD$uniqueHitLabel <- seq_len(nrow(blastsumDD))
#
#cat("Rows after keeping top 10 per scaffold:", nrow(blastsumDD), "\n")
#cat("Unique scaffolds remaining:", length(unique(blastsumDD$Hit)), "\n")
#gc()
#################END OF NEW PRUNING STEP 4/2#######################


#########cat("STEP 6.6: keeping top hits on crowded scaffolds\n")#####commented out 4/2#########old filter step to help with prunning###
#blastsumDD <- blastsumDD[order(
#  blastsumDD$Hit,
#  blastsumDD$E.value,
#  -blastsumDD$Bit.score,
#  -blastsumDD$Hit_length,
#  -blastsumDD$Percent_identity
#), ]
#
#blastsumDD <- do.call(rbind, lapply(split(blastsumDD, blastsumDD$Hit), function(df) {
#  if (nrow(df) > 200) {
#    head(df, 50)
#  } else {
#    df
#  }
#}))
#
#rownames(blastsumDD) <- NULL
#blastsumDD$uniqueHitLabel <- seq_len(nrow(blastsumDD))
#
#cat("Rows after top-hit scaffold filter:", nrow(blastsumDD), "\n")
#gc()

######################END OF PRUNING NEW STEP##################################################################

############################################NEW CODE 3-26###################
cat("STEP 7: scaffold-by-scaffold pruning\n")
cat("Number of rows going into pruning:", nrow(blastsumDD), "\n")

cleanByEValue_prune <- function(grobj, overlaps) {
    queryNums <- unique(queryHits(overlaps))
    keep <- GRanges()

    total_queries <- length(queryNums)
    pb <- txtProgressBar(min = 0, max = total_queries, style = 3)

    while (length(queryNums) > 0) {
        q <- queryNums[1]

        subs_q <- subjectHits(overlaps[queryHits(overlaps) == q])
        grobj_subs <- grobj[subs_q]
        keepThisHit <- which.min(grobj_subs$evalue)
        keep <- c(keep, grobj_subs[keepThisHit])
        queryNums <- queryNums[!(queryNums %in% subs_q)]

        setTxtProgressBar(pb, total_queries - length(queryNums))
    }

    close(pb)
    unique(keep)
}

scaffolds <- unique(blastsumDD$Hit)
cat("Number of unique scaffolds:", length(scaffolds), "\n")

keep_labels <- c()

for (scaf in scaffolds) {
    cat("Processing scaffold:", scaf, "\n")

    scaf_df <- blastsumDD[blastsumDD$Hit == scaf, ]
    cat("Rows on scaffold:", nrow(scaf_df), "\n")

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

    cat("Overlap pairs on scaffold:", length(queryHits(overlaps_scaf)), "\n")

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

pullOut <- blastsumDD[blastsumDD$uniqueHitLabel %in% unique(keep_labels), ]

cat("STEP PRUNING DONE\n")
cat("Number kept after pruning:", nrow(pullOut), "\n")
gc()
#######################################END OF 3-26 CODE#####################
########## Find overlaps, maxgap=100,minoverlap=1,drop hits to itself, drop redundant hits #######
##########################################OLD CODE 
#cat("Number of unique scaffolds:", length(unique(blastsumDD$Hit)), "\n")
#overlaps <- findOverlaps(grobj,maxgap=100L,minoverlap=0L,drop.self=F,drop.redundant=F) 
#queryNums <- unique(queryHits(overlaps))

#cat("STEP OVERLAP CHECK DONE\n")

###############OLD CODE END################################
############################################ RUN PRUNING FUNCTION #######################

#################################################################################
# OLD PRUNE FUNC 
#################################################################################
#cleanByEValue_prune <- function(grobj,overlaps) {
#    queryNums <- unique(queryHits(overlaps))
#    keep <- GRanges()

    # Initialize the progress bar    
#    total_queries <- length(queryNums)
#    pb <- txtProgressBar(min = 0, max = total_queries, style = 3)
#   
#    cat("Number of queryNums entering prune:", total_queries, "\n")
#
#  
#
#    while(length(queryNums) > 0){
#        q <- queryNums[1]
#
#        subs_q <- subjectHits(overlaps[queryHits(overlaps)==q]) # don't need to add the query number in because I have the hits to itself kept into overlaps now     
#        grobj_subs <- grobj[subs_q]
#        keepThisHit <- which.min(grobj_subs$evalue)
#        keep <- c(keep,grobj_subs[keepThisHit])
#        queryNums <- queryNums[!(queryNums %in% subs_q)]
#
#        # print(paste("dropped:",subs_q))
#        # print(paste("left to go: ", length(queryNums))) 
#If there are two ## like for the two prit lines they were ocmmented out in the original copy of the code 
#
#        # Update the progress bar
#        setTxtProgressBar(pb, total_queries - length(queryNums))
#    }
#   
#
#    close(pb)
#    return(unique(keep))
#}

#cat("=== Pruning in process === \n")

#keep_prune_eval <- cleanByEValue_prune(grobj,overlaps)


# check it! should only be self:self matches left! 
# check <- findOverlaps(keep_prune_eval,maxgap=100L,minoverlap=0L,drop.self=F,drop.redundant=F)

# write out, commented out 3-26
#pullOut <- blastsumDD[blastsumDD$uniqueHitLabel %in% keep_prune_eval$uniqueHitLabel,] # pull out the original entries 

# print(head(pullOut))#commented out in original script 

############################# WRITE OUT .sum and .bed #####################################
cat(paste0("=== output files === \npath:", BED_BASE_DIR, "\n"))

write.table(pullOut,
            file.path(BED_BASE_DIR, paste0("cleanedByEvalue.LengthFilterOnly.Pruned.",spp,"_control.sum")),row.names = F,quote=F,sep="/")
cat(" - .sum\n")

pullOut$bed_HitStart0 <- pullOut$Hit_Start - 1 # because bed is zero based (don't change end bc it is non inclusive)
pullOut$bed_HitStart0_minus300 <- pullOut$bed_HitStart0 - 300 # if goes negative, set to 0
pullOut$bed_HitStart0_minus300[pullOut$bed_HitStart0_minus300 < 0] <- 0
pullOut$bed_HitEnd_plus300 <- pullOut$Hit_END + 300 # some may be over limit of chromosome -- deal with that when it happens 


#### check to make sure you don't overhang:

#In Unix run:
# samtools faidx Mustela_putorius_furo.MusPutFur1.0.dna.toplevel.fa &


# need to do samtools faidx on genome:
scaffLengths <- read.table(file.path( SPP_PATH, args[2]),header=F)
colnames(scaffLengths) <- c("Scaffold","length")
pullOut_lengths <- merge(pullOut,scaffLengths[,c("Scaffold","length")],by.x="Hit",by.y="Scaffold",all.x=T,all.y=F)

# check for overhangs:
# pullOut_lengths[pullOut_lengths$bed_HitEnd_plus300 > pullOut_lengths$length,]
# these are the bad ones that need to be fixed

pullOut_lengths[pullOut_lengths$bed_HitEnd_plus300 > pullOut_lengths$length,]$bed_HitEnd_plus300 <- pullOut_lengths[pullOut_lengths$bed_HitEnd_plus300 > pullOut_lengths$length,]$Hit_END

# check for overhangs again:
# pullOut_lengths[pullOut_lengths$bed_HitEnd_plus300 > pullOut_lengths$length,]

pullOut <- pullOut_lengths

##### continue on:
# get strand:
# pullOut$bedstrand <- "." # because bed needs + or - 
# pullOut$bedstrand[pullOut$Hit_strand=="-1"] <- "-"
# pullOut$bedstrand[pullOut$Hit_strand=="1"] <- "+"
# Since Hit_strand has already been written over with '-', '+'. No need to change here. But watch for errors.. 
pullOut$bedstrand <- pullOut$Hit_strand

# for the bedtools name want to have it be like Gang Li's script:
# $species\_$chr\_$start\_$end\_$string\ *** NOTE THAT BLAST IS ONE BASED *** so start end are 1 based
pullOut$bedname <- paste(spp,pullOut$Hit,pullOut$Hit_Start,pullOut$Hit_END,pullOut$Hit_strand,sep="_")
bed <- pullOut[,c("Hit","bed_HitStart0_minus300","bed_HitEnd_plus300","bedname","E.value","bedstrand")]


# SAVE EXPORT FILES 
write.table(bed,file.path(BED_BASE_DIR, paste0(spp, "_control_ORs_OG.cleanedByEvalue.LengthFilterOnly.Pruned.scaff.0based.start.minus300.stop.plus300.name.eval.strand.bed")),row.names = F,col.names = F,quote=F,sep="\t")
cat(" - .bed\n")

write.table(pullOut,file.path(BED_BASE_DIR, paste0(spp, "_control.cleanedByEvalue.noLengthFilter.Pruned.all.pullOut.info.txt")),row.names = F,col.names = T,quote=F,sep="\t")
cat(" - .txt\n")

cat("=== END Rscript Filter_ORs.R ===\n")
