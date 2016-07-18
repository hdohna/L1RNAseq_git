# The following script subset all tophat bam files to retain minimum quality

# Source start script
source('/home/hzudohna/L1polymORFgit/Scripts/_Start_L1polymORF_scg4.R')

# Specifify path to folder that contains all the bam files
FolderPath <- '/srv/gsfs0/projects/levinson/hzudohna/RNAseq/'

# Specify minimum quality value
MinQual <- 10
NewSuffix <- paste("qual", MinQual, ".bam", sep = "")

# Path to script file
ScriptFile <- '/home/hzudohna/qsubMinQual'

# Get all bam files
BamFiles <- list.files(FolderPath, pattern = 'filteredReadOrientation.bam', 
                       recursive = T, full.names = T)
BamFiles <- BamFiles[grep('tophat/', BamFiles)]
BamFiles <- BamFiles[-grep(".bam.", BamFiles)]

# Loop through bam files create a script for indexing per bam file and send off
# script
BamFile <- BamFiles[1]
for (BamFile in BamFiles){

  NameSplit <- strsplit(BamFile, "\\/")[[1]]
  ScriptFileN <- paste(ScriptFile, NameSplit[length(NameSplit) - 2], sep = "_")
  
  # Command line file for sorting
  OutFile <- ReplaceSuffix(BamFile, NewSuffix)
  CommandLine <- c("module load samtools", 
                   paste("samtools view -h -q", MinQual,
                       BamFile, "-o", OutFile))
  cat("Running script", ScriptFileN, "\n")
  CreateAndCallqsubScript(file = ScriptFileN, 
                          qsubCommandLines = CommandLine, 
                          scriptName = 'MinQual')
  
}