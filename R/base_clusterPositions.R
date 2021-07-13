.clusterPositions = function(countFile, chrList, outfileLoc="./", chrom.ref, outfile){
    setwd(outfileLoc)
    ## Calculate number of positions with the same counts
    script <- 'quantile_perl_table_tagcount.pl'
    Fn.Path <- system.file( file.path("Perl",script), package="permseq")
    CMD <- paste("perl ", Fn.Path, " ", countFile, " ", chrom.ref , " ", outfile, "_output_Dnasequantile.txt", sep="")
    system(paste('rm -rf ', outfile, "_output_Dnasequantile.txt", sep=''))
    for(i in chrList){
        CMD <- paste(CMD, i, sep=' ')
    }
    system(CMD, intern = TRUE )
    ## Calculate dnaseThres and knots (using quantileFunction)
    DNase_tags <- read.table(paste(outfileLoc, '/', outfile, '_output_Dnasequantile.txt', sep=''))
    out <- DNase_tags[order(DNase_tags[, 1], decreasing=F), ]
    q <- .quantileFunction(DNase_tags, 0.999)
    q <- c(2, out[1: max(which(out[, 1] <= q)), 1])
    
    cut_off <- NULL
    for(i in 1:length(q)){
        cut_off <- paste(cut_off, q[i], sep=' ')  
    }
    
    Dnase_counts <- c(q, .quantileFunction(DNase_tags, 0.9995))
    dnaseThres <- Dnase_counts
    dnaseKnots <- c(.quantileFunction(DNase_tags, 0.9), .quantileFunction(DNase_tags, 0.99), .quantileFunction(DNase_tags, 0.999))
    for(i in chrList){
    # Partition genome 
        prioroutfile <- paste(i, '_', outfile, '_positions_cluster.txt', sep='')
        script <- 'prior_result.pl'
        Fn.Path <- system.file(file.path("Perl", script), package = "permseq")
        CMD <- paste('perl ', Fn.Path, ' ', chrom.ref, ' ./  ', countFile, ' ', prioroutfile, sep=' ')
        cmd <- paste(CMD, i, cut_off, sep=' ')

        system(cmd, intern = TRUE)

    }
    
    r <- list(dnaseKnots, dnaseThres)
    return(r)
}
