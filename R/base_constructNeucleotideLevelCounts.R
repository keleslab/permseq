.constructNeucleotideLevelCounts = function(fileDir, file, fragL, outfileLoc = "./", chrList, capping){
  system( paste("awk -F \"\\t\" \'{close(f);f=$1}{print > \"", dirname(outfileLoc), "/", basename(outfileLoc), paste("/", "/\"f\"", '_', file, "\"}\' ", sep=''), fileDir, file, sep=""))
  setwd(outfileLoc)
  for(chr.i in chrList){
    file.i <- paste(chr.i, '_', file, sep='')
    .constructBins( infile = file.i, fileFormat = "csem", outfileLoc = outfileLoc, fragLen = fragL, binSize = 1 , capping = capping )
  }
}
