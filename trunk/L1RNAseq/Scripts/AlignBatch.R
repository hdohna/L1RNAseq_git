# The following script get all the paths to fastq files for a batch of RNAseq
# data, creates scripts to align them and sends off the scripts

# Source start script
source('/home/hzudohna/L1polymORF/Scripts/_Start_L1polymORF_bix2.R')

# Specify path to batch folder and pattern for data folder
BatchPath <- '/share/diskarray4/MEI/expression/batch2'
DataPattern <- "sample"

# Specify align command and options
AlignCommand <- '/share/diskarray1/Updated_RNA_software/TopHat/tophat-2.0.9.Linux_x86_64/tophat2' 
AlignOptions <- '--library-type fr-firststrand -o' 

# Specify reference path
RefPath <- '/home/hzudohna/L1polymORF/Data/L1Catalog200FlankSameStrand_Sat_May_07_15-15-31_2016'
ScriptFile <- "/home/hzudohna/qsubScriptAlignRNA"

# Path to output directory
OutputPath <- '/share/diskarray3/hzudohna/RNAseq/'

# Get all paths to fastq files in batch folder
DataDirs <- list.dirs(BatchPath, recursive = F)
DataDirs <- grep(DataPattern, DataDirs, value = T)

# Loop through directories create a qsub file for alignment and submit it
Dir <- DataDirs[1]
for (Dir in DataDirs[-1]){
  FastqFile1 <- list.files(Dir, pattern = '.trimmed.R1.fastq', full.names = T, recursive = T)
  FastqFile2 <- list.files(Dir, pattern = '.trimmed.R2.fastq', full.names = T, recursive = T)
  DirSplit  <- strsplit(Dir, "/")[[1]]
  LastDirPart <- DirSplit[length(DirSplit)]
  OutDir    <-  paste(OutputPath, LastDirPart, sep = "")
  if (!dir.exists(OutDir)) dir.create(OutDir)
  if (length(list.files(OutDir)) == 0){
    CommandLine <- paste(AlignCommand, AlignOptions, OutDir,RefPath, FastqFile1, FastqFile2)
    cat("Running script", ScriptFile, "\n")
    cat("Saving results to", OutDir, "\n")
    CreateAndCallqsubScript(file = ScriptFile, qsubCommandLines = CommandLine, 
                            scriptName = 'AlignRNAseq')
    
  }
}
