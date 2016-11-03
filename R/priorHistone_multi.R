priorHistone_multi = function(object = NULL, dnaseIndex = 1, outfileLoc="./", outfile = "histoneOnly"){

  # Obtain onformation from "object"
  histoneName <- object[['histoneName']]
  chrList <- object[['chrList']]
  chrom.ref <- object[['chrom.ref']]

  s <- scan(chrom.ref, what='char', sep='\n')
  chr_all <- strsplit(s[3], ' ')[[1]]
  chr_length <- as.numeric(strsplit(s[2], ' ')[[1]])
  
  #Generate multipatterns
  message( "Info: Partitioning genome based on multiple data sets..." )
  for(i in 1:length(chrList)){
    chr <- chrList[i]
    file_input <- object@dnaseHistone[[histoneName[dnaseIndex]]][['posLoc_bychr']][[chr]]
    # For Histone files, except the one treated as DNase
    for(h in histoneName[-dnaseIndex]){
      file_input <- paste(file_input, ' ', object@dnaseHistone[[h]][['posLoc3_bychr']][[chr]], sep='')
    }
    # Partition genome based on multiple data sets
    script <- 'multi_data_pattern.pl'
    Fn.Path <- system.file( file.path("Perl", script), package="permseq")
    system(paste('perl ', Fn.Path, ' ', outfileLoc, '/', chr, '_', outfile, '_positions_cluster.txt', " ", chr_length[which(chr_all == chr)], ' ', file_input, sep=''))
  }
  link <- vector('list', length(chrList))
  names(link) <- chrList
  for(i in chrList){
    link[[i]] <- paste(outfileLoc, '/', i, '_', outfile, '_positions_cluster.txt', sep='')
  }
  # Summarizing information
  object@dnaseKnots <- object@dnaseHistone[[histoneName[dnaseIndex]]][['dnaseKnots']]
  object@dnaseThres <- object@dnaseHistone[[histoneName[dnaseIndex]]][['dnaseThres']]
  object@posLoc_bychr <- link             
  object@dnaseName <- object@histoneName[dnaseIndex]
  object@histoneName <- histoneName[-dnaseIndex]
  object@dnaseAlign <- object@histoneAlign[[histoneName[dnaseIndex]]]
  object@histoneNum <- object@histoneNum - 1
  object@histoneAlign <- object@histoneAlign[histoneName[-dnaseIndex]]

  chipSAM <- object@chipSAM
  chipFileName <- gsub(".*/(.*).sam", "\\1", chipSAM)
  outfile_chip <- paste(chipFileName, ".sam", sep = "")
  dnaseThres <- object@dnaseThres
  dnaseKnots <- object@dnaseKnots
  
  outfile_chipmean <- paste("chip", 1:length(outfile_chip), "_chipmean", sep="")
  
  posLoc_bychr <- vector('list', length(object@dnaseName))
  names(posLoc_bychr) <- object@dnaseName
  posLoc_bychr[[length(object@dnaseName)]] <- object@posLoc_bychr

  #calculate averaged chip read counts according to object (dnase or histone information) so that we save time for plotting
  if(!file.exists((paste(outfileLoc, '/', names(posLoc_bychr)[1], '_', outfile_chipmean[1], sep='')))){
    .chipMeanCounts(object, posLoc_bychr, chipSAM, outfile_chip, outfileLoc, outfile_chipmean)
  }
  return(object)
}
