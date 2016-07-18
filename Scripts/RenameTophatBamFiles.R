# The following script renames all tophat bam files so that they can be 
# distinguished when put in the same folder

# Specify path to folder that contains all the bam files
FolderPath <- '/srv/gsfs0/projects/levinson/hzudohna/RNAseq/'

# Get all bam files
BamFiles <- list.files(FolderPath, pattern = 'filteredReadOrientation.qual10.sorted.bam', 
                       recursive = T, full.names = T)
BamFiles <- BamFiles[grep('tophat/', BamFiles)]
BaiFiles <- BamFiles[grep(".bam.bai", BamFiles)]
BamFiles <- setdiff(BamFiles, BaiFiles)
NewBamFileNames <- sapply(BamFiles,function(x){
  PathSplit <- strsplit(x, "/")[[1]]
  FolderPath <- paste(PathSplit[-length(PathSplit)], collapse = "/")
  FolderName <-   PathSplit[length(PathSplit) - 2]
  FileName   <-   PathSplit[length(PathSplit)]
  FileName   <- gsub(".bam", "", FileName)
  paste(FolderPath, "/", FileName, "_", FolderName, ".bam", sep = "")
})
NewBaiFileNames <- sapply(BaiFiles,function(x){
  PathSplit <- strsplit(x, "/")[[1]]
  FolderPath <- paste(PathSplit[-length(PathSplit)], collapse = "/")
  FolderName <-   PathSplit[length(PathSplit) - 2]
  FileName   <-   PathSplit[length(PathSplit)]
  FileName   <- gsub(".bam.bai", "", FileName)
  paste(FolderPath, "/", FileName, "_", FolderName, ".bam.bai", sep = "")
})

file.rename(BamFiles, NewBamFileNames)
file.rename(BaiFiles, NewBaiFileNames)