# The following script renames all tophat bam files so that they can be 
# distinguished when put in the same folder

# Source start script
source('/home/hzudohna/L1polymORFgit/Scripts/_Start_L1polymORF_scg4.R')

# Specifify path to folder that contains all the bam files
FolderPath <- '/srv/gsfs0/projects/levinson/hzudohna/RNAseq'

# Path to script file
ScriptFile <- '/home/hzudohna/qsubIndex'

# Get all bam files
BamFiles <- list.files(FolderPath, pattern = 'filteredReadOrientation.qual10.sorted.bam', 
                       recursive = T, full.names = T)
BamFiles <- BamFiles[grep('tophat/', BamFiles)]
if (length(grep(".bam.", BamFiles)) > 0){
  BamFiles <- BamFiles[-grep(".bam.", BamFiles)]
}

# Loop through bam files create a script for indexing per bam file and send off
# script
for (BamFile in BamFiles){
  
  NameSplit <- strsplit(BamFile, "\\/")[[1]]
  ScriptFileN <- paste(ScriptFile, NameSplit[length(NameSplit) - 2], sep = "_")

  CommandLine <- c("module load samtools", paste('samtools index', BamFile))
  cat("Running script", ScriptFileN, "\n")
  CreateAndCallqsubScript(file = ScriptFileN, 
                          qsubCommandLines = CommandLine, 
                          scriptName = 'SamIndex')
  
}