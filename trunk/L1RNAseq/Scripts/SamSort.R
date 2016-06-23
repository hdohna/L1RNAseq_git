# The following script subset all tophat bam files to retain minimum quality

# Source start script
source('/home/hzudohna/L1polymORF/Scripts/_Start_L1polymORF_bix2.R')

# Specifify path to folder that contains all the bam files
FolderPath <- '/share/diskarray3/hzudohna/RNAseq'

# Path to script file
ScriptFile <- '/home/hzudohna/qsubSort'

# Get all bam files
BamFiles <- list.files(FolderPath, pattern = 'qual10.bam', 
                       recursive = T, full.names = T)
if (length(grep(".bam.", BamFiles)) > 0){
  BamFiles <- BamFiles[-grep(".bam.", BamFiles)]
}

# Loop through bam files create a script for indexing per bam file and send off
# script
for (BamFile in BamFiles){

  NameSplit <- strsplit(BamFile, "\\/")[[1]]
  ScriptFileN <- paste(ScriptFile, NameSplit[length(NameSplit) - 1], sep = "_")

    # Command line file for sorting
  OutFile <- ReplaceSuffix(BamFile, "sorted.bam")
  CommandLine <- paste("/home/txw/samtools/samtools-1.2/samtools sort -o", 
                   OutFile, "-T", paste(BamFile, ".tmp", sep = ""), 
                   BamFile)
  cat("Running script", ScriptFile, "\n")
  CreateAndCallqsubScript(file = ScriptFileN, 
                          qsubCommandLines = CommandLine, 
                          scriptName = 'SamSort')
  
}