.clusterPositions_3cluster = function(countFile, chrList, outfileLoc="./", chrom.ref, outfile){
  setwd(outfileLoc)
  # Calculate number of positions with the same counts
  script <- 'quantile_perl_table_tagcount.pl'
  Fn.Path <- system.file( file.path("Perl", script), package = "permseq")
  CMD <- paste("perl ", Fn.Path, " ", countFile, " ", chrom.ref, " ", outfile, "_output_Dnasequantile.txt", sep = "")
  system(paste('rm -rf ', outfile, "_output_Dnasequantile.txt", sep=''))
  for(i in chrList){
    CMD <- paste(CMD, i, sep = ' ')
  }

  system(CMD, intern = TRUE )
  # Three clusters
  DNase_tags <- read.table(paste(outfileLoc, '/', outfile, '_output_Dnasequantile.txt', sep=''))
  out <- DNase_tags[order(DNase_tags[, 1], decreasing=F), ]
  dnaseKnots <- c(.quantileFunction(DNase_tags, 0.9), .quantileFunction(DNase_tags, 0.99))
  dnaseThres <- dnaseKnots
  cut_off <- NULL
  for (i in 1:length(dnaseThres))
    cut_off <- paste(cut_off, dnaseThres[i], sep = " ")
  for(i in chrList){
    # Partition genome
    prioroutfile <- paste(i, '_', outfile, '_positions_3cluster.txt', sep='')
    script <- 'prior_result.pl'
    Fn.Path <- system.file( file.path("Perl", script), package="permseq")
    CMD <- paste('perl ', Fn.Path, ' ', chrom.ref, ' ./  ', countFile, ' ', prioroutfile, sep=' ')
    cmd <- paste(CMD, i, cut_off, sep=' ')
    system(cmd, intern = TRUE )
    print(i)
  }
}
