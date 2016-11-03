.chipMeanCounts = function(prior = NULL, position_cluster, chipSAM, ChipSamFile, OutfileLoc="./", Outfile){
  setwd(OutfileLoc)
  
  if(length(prior[['chipUni']] == 0)){#if the uni_bed file is not available now
                                        # Obtain uni-reads from each ChIP file  
    for(i in 1:length(ChipSamFile)){
      ChipUniFile <- paste(ChipSamFile[i], '.uni.bed', sep='')
      script <- 'find_uni_bed.pl'
      Fn.Path <- system.file(file.path("Perl", script), package="permseq")   
      CMD <- paste('perl  ', Fn.Path, ' ', chipSAM, ' ', OutfileLoc, '/', ChipUniFile, sep='')
      system(CMD, intern=TRUE)
    }
  }
  # Obtain ChIP mean counts for each data
  for(j in 1:length(position_cluster)){
    CMD <- 'cat '
    for(l in 1:length(position_cluster[[j]])){
      CMD <- paste(CMD, position_cluster[[j]][[l]], sep=' ')
    }
      CMD <- paste(CMD, '>cluster_temp.txt', sep='')
    system(CMD)
    # Calculate averaged ChIP counts at different read counts
    script <- 'mean_chip_update.pl'
    Fn.Path <- system.file(file.path("Perl", script), package="permseq")
    
    for(i in 1:length(ChipSamFile)){
      ChipUniFile <- paste(ChipSamFile[i], '.uni.bed', sep='')
      system(paste('rm -rf ', OutfileLoc, '/', names(position_cluster)[j], '_', Outfile[i], '.txt', sep='') )
      CMD <- paste('perl  ', Fn.Path, ' ', OutfileLoc, '/', ChipUniFile, ' ', 'cluster_temp.txt', ' ', OutfileLoc, '/', names(position_cluster)[j], '_', Outfile[i], sep='')   
      system(CMD,intern=TRUE)
    }
    system('rm -rf cluster_temp.txt')
  }
  
}
