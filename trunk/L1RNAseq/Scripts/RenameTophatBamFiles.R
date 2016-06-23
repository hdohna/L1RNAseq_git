# The following script renames all tophat bam files so that they can be 
# distinguished when put in the same folder

# Specifify path to folder that contains all the bam files
FolderPath <- '/share/diskarray3/hzudohna/RNAseq'

# Get all bam files
BamFiles <- list.files(FolderPath, pattern = 'accepted_hits.bam', 
                       recursive = T, full.names = T)
NewBamFileNames <- sapply(BamFiles,function(x){
  PathSplit <- strsplit(x, "/")[[1]]
  FolderPath <- paste(PathSplit[-length(PathSplit)], collapse = "/")
  FolderName <-   PathSplit[length(PathSplit) - 1]
  FileName   <-   PathSplit[length(PathSplit)]
  FileName   <- gsub(".bam", "", FileName)
  paste(FolderPath, "/", FileName, "_", FolderName, ".bam", sep = "")
})

file.rename(BamFiles, NewBamFileNames)