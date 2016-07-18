# The following script loops through bam files and filters each bam file to get
# only first reads mapped to the reverse stand and second reads mapped to the 
# forward strand. It produces 2 output bam files (one for first reads and one 
# for second reads) per input bam file. These files need to be merged by 
# running MergeFilteredL1Reads.R'

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
InputFileSuffix <- 'accepted_hits.bam'

# Path to script file
ScriptFile <- '/home/hzudohna/qsubFilterR1R2'

# Get all bam files to be filtered
BamFiles <- list.files(BamFolder, pattern = InputFileSuffix, recursive = T, 
                                  full.names = T)
if (length(grep(".bam.", BamFiles)) > 0){
  BamFiles <- BamFiles[-grep(".bam.", BamFiles)]
}

# Loop through bam files and filter them
BamFile <- BamFiles[1]
for (BamFile in BamFiles){
  cat("Filtering bam file", BamFile, "...\n")
  NameSplit <- strsplit(BamFile, "\\/")[[1]]
  ScriptFileN1 <- paste(ScriptFile, NameSplit[length(NameSplit) - 2], "Read1", 
                        sep = "_")
  ScriptFileN2 <- paste(ScriptFile, NameSplit[length(NameSplit) - 2], "Read2", 
                        sep = "_")
  OutFile1 <- gsub("\\.bam", "filteredR1.bam", BamFile)
  OutFile2 <- gsub("\\.bam", "filteredR2.bam", BamFile)
  CommandLine1 <- c("module load samtools", 
                   paste('samtools view -h -f 0x40', BamFile, 
                         '| samtools view -h -f 0x10 -b - >', OutFile1))
  CommandLine2 <- c("module load samtools", 
                    paste('samtools view -h -f 0x80', BamFile, 
                          '| samtools view -h -f 0x20 -b - >', OutFile2))
  cat("Running script", ScriptFileN1, "\n")
  CreateAndCallqsubScript(file = ScriptFileN1, 
                          qsubCommandLines = CommandLine1, 
                          scriptName = 'FilterRead1')
  cat("Running script", ScriptFileN2, "\n")
  CreateAndCallqsubScript(file = ScriptFileN2, 
                          qsubCommandLines = CommandLine2, 
                          scriptName = 'FilterRead2')
  
}

