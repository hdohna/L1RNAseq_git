# The following script get all the paths to fastq files for a batch of RNAseq
# data, creates scripts to align them and sends off the scripts

# Source start script
source('/home/hzudohna/L1polymORFgit/Scripts/_Start_L1polymORF_scg4.R')

# Specify path to batch folder and pattern for data folder
BatchPath <- '/srv/gsfs0/projects/levinson/hzudohna/RNAseq/'
DataPattern <- "sample"

# Specify align command and options
AlignCommand <- c('module load tophat', 'module load bowtie2', 'tophat2')
AlignOptions <- '--library-type fr-firststrand --report-secondary-alignments --max-multihits 1 -o'
# AlignCommand <- c('module load bwa', 'bwa mem')
# AlignOptions <- '-t6' 
AlignTool <- strsplit(AlignCommand[1], 'module load ')[[1]][2]

# Specify reference path
RefPath   <- '/srv/gsfs0/projects/levinson/hzudohna/RefSeqData/L1Catalog200FlankSameStrand_Sat_May_07_15-15-31_2016'
# RefPath   <- '/srv/gsfs0/projects/levinson/hzudohna/RefSeqData/L1Catalog200FlankSameStrand_Sat_May_07_15-15-31_2016.fas'
ScriptDir <- "/home/hzudohna/"

# Path to output directory
OutputPath <- '/srv/gsfs0/projects/levinson/hzudohna/RNAseq/'

# Get all paths to fastq files in batch folder
DataDirs <- list.dirs(BatchPath, recursive = F)
DataDirs <- grep(DataPattern, DataDirs, value = T)

# Loop through directories create a qsub file for alignment and submit it
Dir <- DataDirs[1]
for (Dir in DataDirs){
  FastqFile1 <- list.files(Dir, pattern = '.trimmed.R1.fastq', full.names = T, recursive = T)
  FastqFile2 <- list.files(Dir, pattern = '.trimmed.R2.fastq', full.names = T, recursive = T)
  DirSplit  <- strsplit(Dir, "/")[[1]]
  LastDirPart <- DirSplit[length(DirSplit)]
  OutDir    <-  paste(OutputPath, LastDirPart, sep = "")
  OutDir    <-  paste(OutDir, AlignTool, sep = "/")
  if (!dir.exists(OutDir)) dir.create(OutDir)
  ScriptFile <- paste(ScriptDir, "qsubAlign_", LastDirPart, AlignTool, sep = "")
  if (!dir.exists(OutDir)) dir.create(OutDir)
#  if (length(list.files(OutDir)) == 0){
    CommandLine <- c(AlignCommand[1:2], paste(AlignCommand[3], AlignOptions,
                                            OutDir, RefPath, FastqFile1, FastqFile2))
    # CommandLine <- c(AlignCommand[1],
    #                  paste(AlignCommand[2],RefPath, FastqFile1, FastqFile2,
    #                      ">", paste(OutDir, "/", LastDirPart, "_alnBWA.sam", sep = "")))
    cat("Running script", ScriptFile, "\n")
    cat("Saving results to", OutDir, "\n")
    CreateAndCallqsubScript(file = ScriptFile, 
       qsubHeaderLines = c('#! /bin/sh', '#', '#$ -N TEST', '#', '#$ -cwd', '#',
                           '#$ -l h_rt=48:00:00', '#', '#$ -j y', '#',
                           '#$ -P large_mem', '#', #'#$ -pe shm 6', '#', 
                            '#$ -S /bin/bash', '#', ''),
                qsubCommandLines = CommandLine, 
                    scriptName = paste('AlignRNA', AlignTool, sep = "_"))
    
#  }
}
