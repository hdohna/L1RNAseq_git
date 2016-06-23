# The following script subset all tophat bam files to retain minimum quality

# Source start script
source('/home/hzudohna/L1polymORF/Scripts/_Start_L1polymORF_bix2.R')

# Specifify path to folder that contains all the bam files
FolderPath <- '/share/diskarray3/hzudohna/RNAseq'

# Specify minimum quality value
MinQual <- 10
NewSuffix <- paste("qual", MinQual, ".bam", sep = "")

# Path to script file
ScriptFile <- '/home/hzudohna/qsubMinQual'

# Get all bam files
BamFiles <- list.files(FolderPath, pattern = 'accepted_hits', 
                       recursive = T, full.names = T)
BamFiles <- BamFiles[-grep(".bam.", BamFiles)]

# Loop through bam files create a script for indexing per bam file and send off
# script
for (BamFile in BamFiles){

  NameSplit <- strsplit(BamFile, "\\/")[[1]]
  ScriptFileN <- paste(ScriptFile, NameSplit[length(NameSplit) - 1], sep = "_")
  
  # Command line file for sorting
  OutFile <- ReplaceSuffix(BamFile, NewSuffix)
  CommandLine <- paste("/home/txw/samtools/samtools-1.2/samtools view -h -q", MinQual,
                       BamFile, "-o", OutFile)
  cat("Running script", ScriptFile, "\n")
  CreateAndCallqsubScript(file = ScriptFileN, 
                          qsubCommandLines = CommandLine, 
                          scriptName = 'MinQual')
  
}