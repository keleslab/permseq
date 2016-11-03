.clusterPositions_histoneonly = function(countFile, chrList, outfileLoc="./", chrom.ref, outfile){
  
  setwd(outfileLoc)
  script <- 'quantile_perl_table_tagcount.pl'
  Fn.Path <- system.file( file.path("Perl",script), package="permseq")
  CMD <- paste("perl ", Fn.Path, " ", countFile, " ", chrom.ref, " ", outfile, "_output_Histonequantile.txt", sep="")
  system(paste('rm -rf ', outfile, "_output_Histonequantile.txt", sep=''))
  for(i in chrList){
    CMD <- paste(CMD, i, sep=' ')
  }
  system(CMD, intern = TRUE)

  Histone_tags <- read.table(paste(outfileLoc, '/', outfile, '_output_Histonequantile.txt', sep=''))
  out <- Histone_tags[order(Histone_tags[, 1], decreasing=F), ]
  q <- .quantileFunction(Histone_tags, 0.999)
  q <- c(2, out[1: max(which(out[, 1] <= q)), 1])
  cut_off <- NULL
  for(i in 1:length(q)){
    cut_off <- paste(cut_off, q[i], sep=' ')    
  }

  Histone_counts <- c(q, .quantileFunction(Histone_tags, 0.9995))
  histoneThres <- Histone_counts
  histoneKnots <- c(.quantileFunction(Histone_tags, 0.9), .quantileFunction(Histone_tags, 0.99), .quantileFunction(Histone_tags, 0.999))
  for(i in chrList){
    # Partition genome (Histone)
    prioroutfile <- paste(i, '_', outfile, '_positions_cluster.txt', sep='')
    script <- 'prior_result.pl'
    Fn.Path <- system.file(file.path("Perl", script), package="permseq")
    CMD <- paste('perl ', Fn.Path, ' ', chrom.ref, ' ./  ', countFile, ' ', prioroutfile, sep=' ')
    cmd <- paste(CMD, i, cut_off, sep=' ')
    system(cmd, intern = TRUE )
    #print(i)
  }
  r <- list(histoneKnots, histoneThres)
  # Three clusters
  histoneKnots <- c(.quantileFunction(Histone_tags, 0.9), .quantileFunction(Histone_tags, 0.99))
  histoneThres <- histoneKnots
  cut_off <- NULL
  for (i in 1:length(histoneThres))
    cut_off <- paste(cut_off, histoneThres[i], sep = " ")

  for(i in chrList){
    # Partition genome (Trinary Histone)
    prioroutfile <- paste(i, '_', outfile, '_positions_3cluster.txt', sep='')
    script <- 'prior_result.pl'
    Fn.Path <- system.file(file.path("Perl", script), package="permseq")
    CMD <- paste('perl ', Fn.Path, ' ', chrom.ref, ' ./  ', countFile, ' ', prioroutfile, i,cut_off, sep=' ')
   
    system(CMD, intern = TRUE )
    #print(i)
  }
  return(r)
}
