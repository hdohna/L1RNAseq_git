# The following script renames all tophat bam files so that they can be 
# distinguished when put in the same folder

# Source start script
source('/home/hzudohna/L1polymORFgit/Scripts/_Start_L1polymORF_scg4.R')

# Specifify path to folder that contains all the bam files
FolderPath <- '/srv/gsfs0/projects/levinson/hzudohna/RNAseq/'

# Path to script file
ScriptFile <- '/home/hzudohna/qsubSam2Bam'

# Get all bam files
SamFiles <- list.files(FolderPath, pattern = 'alnBWA.sam', 
                       recursive = T, full.names = T)

# Loop through bam files create a script for indexing per bam file and send off
# script
for (SamFile in SamFiles){
  
  NameSplit <- strsplit(SamFile, "\\/")[[1]]
  ScriptFileN <- paste(ScriptFile, NameSplit[length(NameSplit) - 2], sep = "_")
  BamFile <- gsub("\\.sam", ".bam", SamFile)
  CommandLine <- c("module load samtools", 
                   paste('samtools view', SamFile, '-b -h -o', BamFile))
  cat("Running script", ScriptFileN, "\n")
  CreateAndCallqsubScript(file = ScriptFileN, 
                          qsubCommandLines = CommandLine, 
                          scriptName = 'Sam2Bam')
  
}