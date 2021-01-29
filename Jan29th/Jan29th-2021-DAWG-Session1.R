# DAWG January 29th, 2021: Presented by Mara Cloutier

# Processing/filtering 16S rRNA Illumina amplicon dataset

#1. Install Packages####
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("dada2", "Biostrings", "ShortRead")
         
install.packages("pacman")
pacman::p_load(dada2, Biostrings, ShortRead) 
                     
#2. Set Working directory and Sort Files ####
setwd("~/where your downloaded fastq files are/")
path <- setwd("~/where your downloaded fastq files are/")

fnFs <- sort(list.files(path, pattern = "_R1_001.fastq.gz")
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq.gz")

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`,1)
sample.names

#3. Check for primers ####
FWD <- "GTGYCAGCMGCCGCGGTAA"
REV <- "GGACTACNVGGGTWTCTAAT"

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

# we do not need to run the following three lines of code because I already removed sequences with ambiguous base calls
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs[[1]]))

# now we have to remove the primers using cutadapt
cutadapt <- "/change to wherever cutadapt is located/cutadapt" 
system2(cutadapt, args = "--version")

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC) 
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}
# sanity check!
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
       FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
       REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
       REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

# Break! # Road Block

# 4. Check Complexity and Quality ####
plot(seqComplexity(getSequences(fnFs.cut[1])))
plot(seqComplexity(getSequences(fnRs.cut[1])))

plotQualityProfile(fnFs.cut[1:4])
plotQualityProfile(fnRs.cut[1:4])


#5. Filter/Trim Reads ####
filt_path <- file.path(path, "filtered")
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs.cut, filtFs, fnRs.cut, filtRs, truncLen = c(220,220),
                     maxN=0, maxEE=c(2,2), truncQ = 2, compress = TRUE)

out

#6. Model Error Rates ####
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

#7. Dereplicate Sequences ####
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)

#8. Denoise Sequences ####
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errF, multithread = TRUE)

dadaFs[1]
dadaRs[1]

#9. Merge Forward/Reverse Reads ####
merged <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)

head(mergers[[1]])

seqtab <- makeSequenceTable(merged)

dim(seqtab)

table(nchar(getSequences(seqtab)))

seqtab1 <- seqtab[,nchar(colnames(seqtab)) %in% seq(X, X)]

# 10. Remove chimeric sequences ####

seqtab.nochim <- removeBimeraDenovo(seqtab1, method = "consensus", multithread = TRUE, verbose = TRUE)

dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab1)

# 11. How many Seqs were lost at each step? ####
getN <- function(x) sum(getUniques(x))

track <- cbind(out, sapply(dadaFs, getN), sapply(merged getN), rowSums(seqtab1), rowSums(seqtab.nochim))

colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")

rownames(track) <- sample.names

track

write.table(track, "Processed-16S-Sequencing-Counts.txt")

#12. Assign taxonomy ####
set.seed(128)
taxa <- assignTaxonomy(seqtab.nochim, "~/wherever your silva files are/silva_nr99_v138_wSpecies_train_set.fa.gz", multithread = TRUE)
