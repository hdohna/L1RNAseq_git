# The following script renames all tophat bam files so that they can be 
# distinguished when put in the same folder

# Source start script
source('/home/hzudohna/L1polymORF/Scripts/_Start_L1polymORF_bix2.R')

# Specifify path to folder that contains all the bam files
FolderPath <- '/share/diskarray3/hzudohna/RNAseq'

# Path to script file
ScriptFile <- '/home/hzudohna/qsubIndex'

# Get all bam files
BamFiles <- list.files(FolderPath, pattern = '.qual10..sorted.bam', 
                       recursive = T, full.names = T)
if (length(grep(".bam.", BamFiles)) > 0){
  BamFiles <- BamFiles[-grep(".bam.", BamFiles)]
}

# Loop through bam files create a script for indexing per bam file and send off
# script
for (BamFile in BamFiles){
  
  NameSplit <- strsplit(BamFile, "\\/")[[1]]
  ScriptFileN <- paste(ScriptFile, NameSplit[length(NameSplit) - 1], sep = "_")

  CommandLine <- paste('/home/txw/samtools/samtools-1.2/samtools index', BamFile)
  cat("Running script", ScriptFileN, "\n")
  CreateAndCallqsubScript(file = ScriptFileN, 
                          qsubCommandLines = CommandLine, 
                          scriptName = 'SamIndex')
  
}