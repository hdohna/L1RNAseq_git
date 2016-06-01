# The following script reads a bam file of reads aligned to a catalogue of
# full-length L1 and determines which L1 is present

#############################################
#                                           #
#    Source packages and set parameters     #
#                                           #
#############################################

# Source start script
source('D:/L1polymORF/Scripts/_Start_L1polymORF.r')

# Load packages
library(BSgenome.Hsapiens.UCSC.hg38)
library(csaw)
library(ShortRead)
library(Rsamtools)
library(seqinr)
library(ape)

# Specify file paths
BamFolder        <- 'D:/L1RNAseq/Data/'
CatalogueFile  <- "D:/L1polymORF/Data/L1Catalogue_Updated_Sat_May_07_15-15-31_2016.csv"
CatalogSeqFile <- "D:/L1polymORF/Data/L1CatalogueWithFlank_Sat_May_07_15-15-31_2016.fas"
AlignListFile  <- "D:/L1polymORF/Data/L1CatalogueWithFlank_Sat_May_07_15-15-31_2016_L1Locations.RData"

############################
#                          #
#        Read Data         #
#                          #
############################

# Load file with  
load(file = AlignListFile)
if(!all(colnames(L1StartEnd) == names(L1withFlank))){
  stop("Names L1StartEnd and L1withFlank don't match!")
}
# Get accession numbers of sequences that contain L1 
blnSeqWithL1     <- L1StartEnd["Width", ] >= 5000
min(L1StartEnd["Width", ])

# Read catalog table
L1Catalog <- read.csv(CatalogueFile, as.is = T)
AccMatch  <- match(colnames(L1StartEnd), L1Catalog$Accession)
L1Catalog <- L1Catalog[AccMatch, ] 

# Get file paths to all bam files
BamFiles <- list.files(BamFolder, pattern = ".bam", full.names = T)
BamFiles <- BamFiles[-grep(".bam.", BamFiles)]

############################
#                          #
#   Calculate coverage     #
#                          #
############################

# Extract experiment names from file names
ExperimentNames <- sapply(BamFiles, function(x){
  Split1 <- strsplit(x, "/")[[1]]
  strsplit(Split1[length(Split1)], "_")[[1]][1]
})

# Process different bam files
CoverObjectList <- lapply(BamFiles, function(x){
  cat(" *****    Processing", x, "****\n\n")
  CalculateCoverMatL1Catalog(x, L1Catalog,L1withFlank, L1StartEnd)
})
names(CoverObjectList) <- ExperimentNames

# Get total coverage per L1
CoverByExperiment <- sapply(CoverObjectList, function(x)rowSums(x$CoverMatL1))
CoverSum <- rowSums(CoverByExperiment)
plot(L1Catalog$Activity, CoverSum)

############################
#                          #
#     Plot coverage        #
#                          #
############################

# Plot the first 9 elements with highest coverage sums
idxHighCover <- order(CoverSum, decreasing = T)[1:9]

Cols <- rainbow(length(CoverObjectList))
par(mfrow = c(3, 3), oma = c(2, 2, 0, 2))
for (i in idxHighCover){
  plot(CoverObjectList[[1]]$CoverMatL1[i,], xlab = "", 
       ylab = "",
       main = rownames(CoverObjectList[[1]]$CoverMatL1)[i],
       type = "l", col = Cols[1])
  for (j in 2:length(CoverObjectList)){
    lines(CoverObjectList[[j]]$CoverMatL1[i,], col = Cols[j])
  }
}
legend("topleft", legend = ExperimentNames, col = Cols, 
       lty = rep(1, length(ExperimentNames)),
       y.intersp = 0.15, cex = 0.6, bty = "n")
mtext("Position on L1", side = 1, outer = T)
mtext("Coverage", side = 2, outer = T)
CreateDisplayPdf("D:/L1RNAseq/Figures/iPScoverPerL1.pdf")

# Determine peaks and positions for first experiment
ExperimentNames[1]
Exp1 <- 1
CMat <- CoverObjectList[[Exp1]]$CoverMatL1
PeakPos <- lapply(1:nrow(CMat), function(i){
  Peaks <- slice(CMat[i,], lower = 2)
  viewWhichMaxs(Peaks)
})
NrPeaks <- sapply(PeakPos, length)

# Plot different coverages against each other
i    <- 1
Exp1 <- 1
Exp2 <- 3
par(mfrow = c(3, 3), oma = c(2, 2, 0, 2))
for (i in idxHighCover){
  CM1 <- CoverObjectList[[Exp1]]$CoverMatL1
  CM2 <- CoverObjectList[[Exp2]]$CoverMatL1
  plot(CM1[i,],CM2[i,],
       xlab = ExperimentNames[Exp1], ylab = ExperimentNames[Exp2],
       main = rownames(CM1)[i])
  points(CM1[i, PeakPos[[i]]], CM2[i, PeakPos[[i]]], col = "red", pch = 16)
  lines(c(0, 20000), c(0, 20000))
  
}
CreateDisplayPdf("D:/L1RNAseq/Figures/iPScoverScatter.pdf")

#################################
#                               #
#   Calculate 'mapability'      #
#                               #
#################################

# Read alignment of L1s
L1Alignment <- read.dna("D:/L1polymORF/Data/L1Catalogue_aligned.fas", 
                        format = "fasta", as.character = T)

# Get accession and allele numbers and match to other data
NameParts <- t(sapply(names(L1Alignment), function(x) strsplit(x, "_")[[1]]))
colnames(NameParts) <- c("Accession", "Allele")
idxAllele1 <- which(NameParts[,"Allele"] == "1")
AccMatch2  <- match(L1Catalog$Accession, NameParts[idxAllele1,"Accession"])

# Create an alignment matrix and check that accession names match
L1Alignment <- t(sapply(idxAllele1[AccMatch2], function(i) L1Alignment[[i]]))
NameParts2 <- t(sapply(rownames(L1Alignment), function(x) strsplit(x, "_")[[1]]))
colnames(NameParts2) <- c("Accession", "Allele")
all(L1Catalog$Accession == NameParts2[ ,"Accession"])

# Loop through all L1s and create a vector for each L1 that indicates for each
# bp whether it contains a nucleotide unique to this L1
rowIndices <- 1:nrow(L1Alignment)
x <- 1
L1UniqueNucList <- lapply(rowIndices, function(x){
  Seq <- L1Alignment[x,]
  idxNoGap <- Seq != "-"
  Seq <- Seq[idxNoGap]
  DiffMat <- t(sapply(rowIndices[rowIndices != x], 
                    function(y) Seq != L1Alignment[y, idxNoGap]))
  1 * apply(DiffMat, 2, all)
})
plot(L1UniqueNucList[[1]])
plot(L1UniqueNucList[[5]])

NrUnique <- sapply(L1UniqueNucList, sum)
hist(sapply(L1UniqueNucList, sum), breaks = 0:250, xlim = c(0, 50))

plot(NrUnique, CoverSum)
Lmfit <- lm(CoverSum ~ NrUnique)
plot(L1Catalog$Activity, Lmfit$residuals)

# Calculate minimum distance to closest sequence per 100 bp window
StartVals <- seq(1, 6150, 10)
MidPoints <- StartVals + 50
MindistPerWindow    <- matrix(nrow = nrow(L1Alignment), ncol = length(StartVals))
idxMindistPerWindow <- MindistPerWindow
NrMindistPerWindow <- MindistPerWindow
for(i in 1:length(StartVals)) {
    x <- StartVals[i]
    W <- x:(x+100)
    DistMat <- dist.dna(as.DNAbin(L1Alignment[,W]), model = "raw", as.matrix = T, 
                        pairwise.deletion = T)
    diag(DistMat) <- NA
    #    cat("Dimensions of dist:", dim(DistMat), "\n")
    if(! all(colnames(DistMat) == rownames(L1Alignment))){
      stop("Distance matrix not consistent with alignment")
    }
    for (j in 1:ncol(DistMat)){
      MinDist <- min(DistMat[,j], na.rm = T)
      idxMin  <- which(DistMat[,j] == MinDist)
#      cat("Index length is", length(idxMin), "\n")
      idxMindistPerWindow[j,i] <- idxMin[1]
      NrMindistPerWindow       <- length(idxMin)
      MindistPerWindow[j,i]    <- MinDist
    }
}

# Determine for each sequence a map from alignment position to sequence 
# position. The matrix generated below gives for each L1 (row) an index 
# that gives for the respective MidPoints position (column) its position
# in the sequence without indels
StarVals2SeqMap <- sapply(1:nrow(L1Alignment) function(i){
  SeqPos <- which(L1Alignment[i,] != "-")
  which(SeqPos %in% MidPoints)
  
})
# Plot coverage with
par(mfrow = c(3, 3), oma = c(2, 2, 0, 4))
for (i in idxHighCover){
  plot(CoverObjectList[[1]]$CoverMatL1[i,], xlab = "", 
       ylab = "",
       main = rownames(CoverObjectList[[1]]$CoverMatL1)[i],
       type = "l", col = Cols[1])
  par(new = TRUE)
  plot(MidPoints, MindistPerWindow[i,], type = "l", axes = FALSE, bty = "n", 
       xlab = "", ylab = "")
  axis(side=4)
}
mtext("Position on L1", side = 1, outer = T)
mtext("Coverage", side = 2, outer = T)
mtext("Proportional difference to closest sequence", side = 4, outer = T,
      line = 2)
CreateDisplayPdf("D:/L1RNAseq/Figures/iPScoverPerL1withMapability.pdf")

# Determine closely related L1 in 4000-4500 bp window
idxWindow <- MidPoints >= 4000 & MidPoints <= 4500
#####################################
#                                   #
#   Plot coverage & mappability     #
#                                   #
#####################################

# 
# Plot different coverages against each other
Row  <- 1
Exp1 <- 1
Exp2 <- 2
Cover1 <- CoverObjectList[[Exp1]]$CoverMatL1[Row,]
Cover2 <- CoverObjectList[[Exp2]]$CoverMatL1[Row,]
plot(Cover1, Cover2, xlab = ExperimentNames[Exp1], ylab = ExperimentNames[Exp2],
     main = rownames(CoverObjectList[[Exp2]]$CoverMatL1)[Row])
lines(c(0, 200), c(0, 200))
blnUnique <- L1UniqueNucList[[Row]] == 1
sum(blnUnique)
points(Cover1[blnUnique], Cover2[blnUnique], pch = 15, col = "red")

