##############################################
#
# General description:
#
#   The following script reads a repeat table downloaded from the genome
#   browser repeatMasker track (http://genome.ucsc.edu/cgi-bin/hgTables)
#   and gets all old L1s and generates reads from the sequences

# Input:
#
#    /srv/gsfs0/projects/levinson/hzudohna/RefSeqData/repeatsHg38: table with all repeats
#   

# Output:
#   
#    SimulatedOldL1.fastq

##############################################

########################################
#                                      #
#  Source packages and set parameters  #
#                                      #
########################################

# Source start script
#source('D:/L1polymORF/Scripts/_Start_L1polymORF.r')
source('/home/hzudohna/RScripts/_Start_SCG4.R')

# Load packages
library(BSgenome.Hsapiens.UCSC.hg38)

# Files and folders
RepeatFile <- "D:/L1polymORF/Data/repeatsHg38"
RepeatFile <- "/srv/gsfs0/projects/levinson/hzudohna/RefSeqData/repeatsHg38"
OutFastQfilePath1 <- '/srv/gsfs0/projects/levinson/hzudohna/RNAseq/SimulatedOldL1_R1.fastq'
OutFastQfilePath2 <- '/srv/gsfs0/projects/levinson/hzudohna/RNAseq/SimulatedOldL1_R2.fastq'

# Read parameters
MeanCover         <- 40
ReadLength        <- 150
InsertSize        <- 300
InsertSD          <- 50
FragmentSize      <- InsertSize + 2 * ReadLength
ReadsSampledPerL1 <- round(MeanCover * 6000 / ReadLength / 2)

#######################################
#                                     #
#   Read data and get L1 sequences    #
#                                     #
#######################################

cat("Read and process repeatMasker table \n")

# Read repeat table
RepeatTable <- read.delim(RepeatFile)
RepeatTable <- RepeatTable[nchar(as.character(RepeatTable$genoName)) <= 5, ]

# Subset repeat table to get only L1
RepeatTable <- RepeatTable[RepeatTable$repFamily == "L1", ]

# Subset to get only full-length elements
blnFullLength <- abs(RepeatTable$genoEnd - RepeatTable$genoStart) >= 6000
RepeatTable <- RepeatTable[blnFullLength, ]

# Make some corrections to create a proper GRanges object with L1 Seqences
L1IRanges <- IRanges(start = RepeatTable$genoStart, end = RepeatTable$genoEnd)
L1GRanges <- GRanges(seqnames = RepeatTable$genoName, ranges = L1IRanges,
                     strand = RepeatTable$strand)

# Get all L1 sequences  
cat("Get sequences for L1 element \n")
L1Seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38, L1GRanges)

#######################################
#                                     #
#   Read data and get L1 sequences    #
#                                     #
#######################################

# Initialize connections for fastq files for read 1 and 2 and counting 
# variables
OF1 <- file(OutFastQfilePath1, "w")
OF2 <- file(OutFastQfilePath2, "w")
ReadCounter <- 0
QString     <- paste(rep('~', ReadLength), collapse = "")

# Loop over L1s and sample random reads for each L1
i <- 1
for (i in 1: length(L1Seq)){
  L1 <- L1Seq[i]
  if (RepeatTable$strand[i] == "-") L1 <- reverseComplement(L1)
  
  # Sample random paired reads
  R2Start <- sample(width(L1), ReadsSampledPerL1, replace = T)
  R1End   <- R2Start + round(rnorm(ReadsSampledPerL1, FragmentSize, InsertSD))
  blnInL1 <- R1End <= width(L1)
  R2Start <- R2Start[blnInL1]
  R1End   <- R1End[blnInL1]
  L1Rep   <- rep(L1, sum(blnInL1))
  R2Seq   <- subseq(L1Rep, start = R2Start, width = ReadLength)
  R1Seq   <- subseq(L1Rep, end = R1End, width = ReadLength)
  R1Seq   <- reverseComplement(R1Seq)
  SeqLines1 <- as.character(R1Seq)
  SeqLines2 <- as.character(R2Seq)
  NameLines <- paste("@", ReadCounter + (1 : length(SeqLines1)), sep = "")
  ReadLines1 <- rep("+", 4 * length(SeqLines1))
  ReadLines2 <- rep("+", 4 * length(SeqLines2))
  idxNames <- seq(1, length(ReadLines1), 4)
  
  # Add names
  ReadLines1[idxNames] <- NameLines
  ReadLines2[idxNames] <- NameLines
  
  # Add sequences
  ReadLines1[idxNames + 1] <- SeqLines1
  ReadLines2[idxNames + 1] <- SeqLines2
  
  # Add read quality
  ReadLines1[idxNames + 3] <- QString
  ReadLines2[idxNames + 3] <- QString
  
  # Write reads into a fasta file
  writeLines(ReadLines1, OF1)
  writeLines(ReadLines2, OF2)
  ReadCounter <- ReadCounter + length(R1Seq)
  cat("Total of", WriteCounter, " reads written out \n")
  
}
close(OF1)
close(OF1)
