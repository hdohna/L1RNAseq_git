# The following script subset all tophat bam files to retain minimum quality

# Source start script
source('/home/hzudohna/L1polymORFgit/Scripts/_Start_L1polymORF_scg4.R')

# Specifify path to folder that contains all the bam files
FolderPath <- '/srv/gsfs0/projects/levinson/hzudohna/RNAseq'


# Path to script file
ScriptFile <- '/home/hzudohna/qsubSort'

# Get all bam files
BamFiles <- list.files(FolderPath, pattern = 'filteredReadOrientation.qual10.bam', 
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

    # Command line file for sorting
  OutFile <- ReplaceSuffix(BamFile, "sorted.bam")
  CommandLine <- c("module load samtools", paste("samtools sort -o", 
                   OutFile, "-T", paste(BamFile, ".tmp", sep = ""), 
                   BamFile))
  cat("Running script", ScriptFileN, "\n")
  CreateAndCallqsubScript(file = ScriptFileN, 
                          qsubCommandLines = CommandLine, 
                          scriptName = 'SamSort')
  
}