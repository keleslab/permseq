.summary = function(outfileLoc, fileName)
{
  con <- file(paste(outfileLoc, "/", fileName, sep = ""), open = "r")
  line <- readLines(con, 1)
  while(length(grep("# reads processed:", line)) == 0){
    
    line <- readLines(con, 1)
  }
  readNum <- as.numeric(gsub("#.*: (.*)", "\\1", line))

  while(length(grep("# reads with at least one reported alignment", line)) == 0){
        line <- readLines(con, 1)
  }
  readAligned <- as.numeric(gsub("#.*: (.*) .*%)", "\\1", line))
  readAlignedRate <- as.numeric(strsplit(strsplit(line, "\\(")[[1]][2], "%)")[[1]][1])

  while(length(grep("# reads that failed to align", line)) == 0){
    line <- readLines(con, 1)
  }
  readFailed <- as.numeric(gsub("#.*: (.*) .*%)", "\\1", line))
  readFailedRate <- as.numeric(strsplit(strsplit(line, "\\(")[[1]][2], "%)")[[1]][1])
    
  while(length(grep("Reported", line)) == 0){
    line <- readLines(con, 1)
  }

  AlignNum <- as.numeric(gsub("Reported (.*) al.*", "\\1", line))
  close(con)
  system(paste("rm -r -f ", outfileLoc, "/", fileName, sep = ""))
  return(list(readNum = readNum,
              readAlignedNum = readAligned,
              readAlignedRate = readAlignedRate,
              readFailedNum = readFailed,
              readFailedRate = readFailedRate,
              AlignNum = AlignNum))
}
