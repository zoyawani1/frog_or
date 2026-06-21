#!/usr/bin/env Rscript

#################################################################################
# WRITE FUNCTION COMMENTS HERE 
#################################################################################

require(GenomicRanges, quietly = TRUE)
# try(library(GenomicRanges), silent=TRUE)

args = commandArgs(trailingOnly=TRUE)

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

# Set arguments to variables 
SPP_PATH=args[1]
setwd(SPP_PATH)
spp=tail(unlist(strsplit(SPP_PATH, '/')), 1) # Grabs the last directory as spp, should be Specie

cat("=== BEGIN Rscript Filter_ORs.R ===\n\n")

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

blast_sum_file <- file.path(BLAST_DIR, paste0(genome_id, "_OR_query.sum"))

if (!file.exists(blast_sum_file)) {

    cat("WARNING: BLAST sum file not found for genome:", genome_id, "\n")

    write(
        genome_id,
        file = "/project/stuckert/users/wani/f2/anuran_genomes/missing_sum_files.txt",
        append = TRUE
    )

    quit(save="no")
}

cat("=== Using BLAST sum file ===\n", blast_sum_file, "\n")

blastsumDD_raw <- read.table(
  blast_sum_file,
  header = TRUE,
  sep = "/",
  quote = ""
)


# do some filtering ahead of time to speed up processing:
# get rids of hits < 250bp in length
# tail(blastsumDD_raw)
# dim(blastsumDD_raw)

blastsumDD_raw$Hit_length <- as.numeric(as.character(blastsumDD_raw$Hit_length))
blastsumDD <- blastsumDD_raw[blastsumDD_raw$Hit_length >= 250,]
blastsumDD <- na.omit(blastsumDD)
blastsumDD$Hit_Start <- as.numeric(as.character(blastsumDD$Hit_Start))
blastsumDD$Hit_END <- as.numeric(as.character(blastsumDD$Hit_END))

# rownames(blastsumDD) <- NULL
blastsumDD$uniqueHitLabel <- seq(1:dim(blastsumDD)[1])

# Converts strands to '+' and '-'
blastsumDD$Hit_strand <- ifelse(blastsumDD$Hit_strand == ' 1', 
                                '+', 
                                ifelse(blastsumDD$Hit_strand == ' -1', 
                                        '-', blastsumDD$Hit_strand))

# dim(blastsumDD)


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

grobj <- GRanges(seqnames=blastsumDD$Hit,
                 ranges=IRanges(start=blastsumDD$Hit_Start,
                                end=blastsumDD$Hit_END),
                 strand=blastsumDD$Hit_strand,
                 bit.score=blastsumDD$Bit.score,
                 evalue=blastsumDD$E.value,
                 query=blastsumDD$Query,
                 hit.length=blastsumDD$Hit_length,
                 uniqueHitLabel=blastsumDD$uniqueHitLabel)


########## Find overlaps, maxgap=100,minoverlap=1,drop hits to itself, drop redundant hits #######
overlaps <- findOverlaps(grobj,maxgap=100L,minoverlap=0L,drop.self=F,drop.redundant=F) 
queryNums <- unique(queryHits(overlaps))


############################################ RUN PRUNING FUNCTION #######################

#################################################################################
# WRITE FUNCTION DESCRIPTION 
#################################################################################
cleanByEValue_prune <- function(grobj,overlaps) {
    queryNums <- unique(queryHits(overlaps))
    keep <- GRanges()

    # Initialize the progress bar    
    total_queries <- length(queryNums)
    pb <- txtProgressBar(min = 0, max = total_queries, style = 3)
  

    while(length(queryNums) > 0){
        q <- queryNums[1]

        subs_q <- subjectHits(overlaps[queryHits(overlaps)==q]) # don't need to add the query number in because I have the hits to itself kept into overlaps now     
        grobj_subs <- grobj[subs_q]
        keepThisHit <- which.min(grobj_subs$evalue)
        keep <- c(keep,grobj_subs[keepThisHit])
        queryNums <- queryNums[!(queryNums %in% subs_q)]

        # print(paste("dropped:",subs_q))
        # print(paste("left to go: ", length(queryNums)))

        # Update the progress bar
        setTxtProgressBar(pb, total_queries - length(queryNums))
    }
    return(unique(keep))
}

cat("=== Pruning in process === \n")

keep_prune_eval <- cleanByEValue_prune(grobj,overlaps)


# check it! should only be self:self matches left! 
# check <- findOverlaps(keep_prune_eval,maxgap=100L,minoverlap=0L,drop.self=F,drop.redundant=F)

# write out
pullOut <- blastsumDD[blastsumDD$uniqueHitLabel %in% keep_prune_eval$uniqueHitLabel,] # pull out the original entries 

# print(head(pullOut))

############################# WRITE OUT .sum and .bed #####################################
cat(paste0("=== output files === \npath:", SPP_PATH,  "\n"))

write.table(pullOut,
            file.path(SPP_PATH, paste0("cleanedByEvalue.LengthFilterOnly.Pruned.",spp,"_control.sum")),row.names = F,quote=F,sep="/")
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
write.table(bed,file.path(SPP_PATH, paste0(spp, "_control_ORs_OG.cleanedByEvalue.LengthFilterOnly.Pruned.scaff.0based.start.minus300.stop.plus300.name.eval.strand.bed")),row.names = F,col.names = F,quote=F,sep="\t")
cat(" - .bed\n")

write.table(pullOut,file.path(SPP_PATH, paste0(spp, "_control.cleanedByEvalue.noLengthFilter.Pruned.all.pullOut.info.txt")),row.names = F,col.names = T,quote=F,sep="\t")
cat(" - .txt\n")

cat("=== END Rscript Filter_ORs.R ===\n")
