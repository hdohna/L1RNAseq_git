# The following script loops through bam files and merges filtered read1 and 2 
# files created by script FilterL1RNAseqByOrientation.R

#############################################
#                                           #
#    Source packages and set parameters     #
#                                           #
#############################################

# Source start script
source('/home/hzudohna/L1polymORFgit/Scripts/_Start_L1polymORF_scg4.R')

# Specify file paths
BamFolder <- '/srv/gsfs0/projects/levinson/hzudohna/RNAseq/'

# Specify file name parts
InputFileSuffix <- 'filteredR1.bam'
PathPart <- 'tophat/'

# Path to script file
ScriptFile <- '/home/hzudohna/qsubMergeR1R2'

# Get all bam files to be filtered
BamFiles <- list.files(BamFolder, pattern = InputFileSuffix, recursive = T, 
                                  full.names = T)
BamFiles <- BamFiles[grep(PathPart, BamFiles)]
if (length(grep(".bam.", BamFiles)) > 0){
  BamFiles <- BamFiles[-grep(".bam.", BamFiles)]
}

# Loop through bam files and filter them
BamFile <- BamFiles[1]
for (BamFile in BamFiles){
  cat("Filtering bam file", BamFile, "...\n")
  NameSplit <- strsplit(BamFile, "\\/")[[1]]
  ScriptFileN <- paste(ScriptFile, NameSplit[length(NameSplit) - 2], "Merge", 
                        sep = "_")
  BamFile2 <- gsub("filteredR1.bam", "filteredR2.bam", BamFile)
  OutFile <- gsub("filteredR1.bam", "filteredReadOrientation.bam", BamFile)
  CommandLine <- c("module load samtools", 
                   paste('samtools merge', OutFile, BamFile, BamFile2))
  cat("Running script", ScriptFileN, "\n")
  CreateAndCallqsubScript(file = ScriptFileN, 
                          qsubCommandLines = CommandLine, 
                          scriptName = 'MergeReads')
}

