priorProcess = function(dnaseFile = NULL, histoneFile = NULL, dnaseName = 'dnase', histoneName = NULL, fragL = 200, AllocThres = 900, chrList = NULL, capping = 0, outfileLoc = "./", outfile = "dnase_histone", bowtieDir = NULL, bowtieIndex = NULL, vBowtie = 2, mBowtie = 99, pBowtie = 8, bwaDir = NULL, bwaIndex = NULL, nBWA = 2, oBWA = 1, tBWA = 8, mBWA = 99, csemDir = NULL, picardDir = NULL, chrom.ref=NULL, saveFiles = TRUE){


  #Refine the input parameters -- NULL NA numeric ""
  #paraList <- c(dnaseFile, histoneFile, dnaseName, histoneName, fragL, AllocThres, chrList, capping, outfileLoc, outfile, bowtieDir, bowtieIndex, vBowtie, mBowtie, pBowtie, bwaDir, bwaIndex, nBWA, oBWA, tBWA, mBWA, csemDir, chrom.ref, saveFiles)
  #for(para in paraList){
  #  para <- .refineParameters(para)
  #}

  
  # Check parameters
  if(length(dnaseFile)>1){
    stop("Only one DNase-seq file is acceptable. DNase file needs to be merged prior running the priorProcess function.")
  }else if(is.null(dnaseFile)){
    stop("Please provide proper path to DNase-seq file.")
  }
  
  if(is.null(histoneName))
    histoneName <- 1:length(histoneFile)
 
  if(fragL<=0)
    stop("fragL (average fragment length) must be a positive value.")

  if(AllocThres<0 | AllocThres>1000)
    stop("AllocThres (allocation score threshold) must be a number between 0 and 1000.")
  
  if(tolower(sub(".*\\.", "", dnaseFile)) == "fastq" & is.null(bowtieIndex) & is.null(bwaIndex))
    stop("To process DNase-seq fastq file, please provide either Bowtie index or BWA index to accomplish the alignment.")

  if(!(tolower(sub(".*\\.", "", dnaseFile)) %in% c("fastq", "sam", "bam", "bed")) )
    stop("DNase-seq file must be in fastq, sam , bam or bed format.")

  if(prod(tolower(sub(".*\\.", "", histoneFile)) %in% c("fastq", "sam", "bam", "bed")) == 0 )
    stop("All Histone ChIP-seq files must be in fastq, sam , bam or bed format.")

  dnaseFormat <- tolower(sub(".*\\.", "", dnaseFile))
  if(dnaseFormat == "fastq" & is.null(bowtieIndex) & is.null(bwaIndex))
    stop("Please choose one alignment method by providing the corresponding index.")
  
  if(vBowtie<0 | vBowtie>3)
    stop("vBowtie (-v option in Bowtie) must be 0, 1, 2 or 3.")
  
  if(mBowtie%%1!=0 | mBowtie<0)
    stop("mBowtie (-m option in Bowtie) must be a positive integer.")
  
  if(pBowtie%%1!=0 | pBowtie<0)
    stop("pBowtie (-p option in Bowtie) must be a positive integer.")


  if(nBWA < 0 | nBWA > 3)
    stop("nBWA (bwa aln -n option in BWA) must be 0, 1, 2 or 3.")

  if(oBWA < 0 | oBWA > 1)
    stop("oBWA (bwa aln -o option in BWA) must be 0 or 1.")
  
  if(mBWA%%1!=0 | mBWA<0)
    stop("mBWA (bwa samse -n option in BWA) must be a positive integer.")
  
  if(tBWA%%1!=0 | tBWA<0)
    stop("tBWA (bwa aln -t option in BWA) must be a positive integer.")
  
  
  
  #build the reference chromosome file using bowtie-inspect
  if(!dir.exists(outfileLoc))
    system(paste('mkdir ',outfileLoc,sep=''))

  #change into the output folder directory
  setwd(outfileLoc)
  if(is.null(chrom.ref)){
    if(!is.null(bowtieIndex)){
      script <- "genRef.pl"
      Fn.Path <- system.file(file.path("Perl", script), package = "permseq")
      bowtieIndexName <-  gsub(".*/(.*)", "\\1", bowtieIndex)
      chrom.ref <- paste(outfileLoc, "/", bowtieIndexName, ".ref", sep = "")
      CMD <- paste(bowtieDir, "/bowtie-inspect -s ", bowtieIndex, " | perl ", Fn.Path, " ", chrom.ref, sep = "")
      system(CMD, intern = TRUE)
    }
  }

  if(!is.null(bowtieIndex)){
    # Sumarize Bowtie Information
    bowtieInfo <- vector('list', 5)
    names(bowtieInfo) <- c("bowtieIndex","bowtieDir","vBowtie","mBowtie","pBowtie")
    bowtieInfo$bowtieIndex <- bowtieIndex
    bowtieInfo$bowtieDir <- bowtieDir
    bowtieInfo$vBowtie <- vBowtie
    bowtieInfo$mBowtie <- mBowtie
    bowtieInfo$pBowtie <- pBowtie

    bwaInfo <- NULL
  }else{
    
    # Summarize BWA Information
    bwaInfo <- vector('list', 6)
    names(bwaInfo) <- c("bwaIndex","bwaDir","nBWA","oBWA","tBWA", "mBWA")
    bwaInfo$bwaIndex <- bwaIndex
    bwaInfo$bwaDir <- bwaDir
    bwaInfo$nBWA <- nBWA
    bwaInfo$oBWA <- oBWA
    bwaInfo$tBWA <- tBWA
    bwaInfo$mBWA <- mBWA

    bowtieInfo <- NULL
  }
  # if only one DNase data
  if(is.null(histoneFile) & !is.null(dnaseFile)){ 
    # Process DNase data
    dnaseObject <- .dnaseProcess(dnaseFile, fragL, AllocThres, chrList, chrom.ref, capping, outfileLoc = paste(outfileLoc, '/', dnaseName, '/', sep=''), paste(dnaseName, '_', outfile, sep=''), bowtieIndex, csemDir, bowtieDir, vBowtie, mBowtie, pBowtie, bwaDir, bwaIndex, nBWA, oBWA, tBWA, mBWA, picardDir, saveFiles)
    chrList <- dnaseObject[['chrList']]
    chrom.ref <- dnaseObject[['chrom.ref']]
    if( dnaseFormat %in% c("sam", "bam", "bed") || is.null(bowtieIndex)){
      dnaseAlign <- list()  ## If DNase file is sam, bam or bed, or align by BWA: No aligment information
    }else{
      dnaseAlign <- .summary(paste(outfileLoc, '/', dnaseName, sep = ""), "priorProcessDNaseBowtie_temp.txt") ##fastq file aligned by Bowtie will have summary information
    }
    return(new("Prior",
               dnaseKnots = dnaseObject[['dnaseKnots']],                 
               dnaseThres = dnaseObject[['dnaseThres']],                 
               posLoc_bychr= dnaseObject[['posLoc_bychr']],                 
               dnaseAlign = dnaseAlign,
               chrList = chrList,
               dataNum = 1,                 
               dnaseName = dnaseName,                 
               fragL = fragL,                 
               bowtieInfo = bowtieInfo,
               bwaInfo = bwaInfo,
               csemDir = csemDir,          
               picardDir = picardDir,
               outfileLoc = outfileLoc,             
               chrom.ref = chrom.ref))
    
     
  }else if(!is.null(histoneFile) & !is.null(dnaseFile)){# If DNase and Histone data
    histoneNum <- length(histoneFile)
    # Process DNase data
    dnaseObject <- .dnaseProcess(dnaseFile, fragL, AllocThres, chrList, chrom.ref, capping, outfileLoc = paste(outfileLoc, '/', dnaseName, '/', sep=''), paste(dnaseName, '_', outfile, sep=''), bowtieIndex, csemDir, bowtieDir, vBowtie, mBowtie, pBowtie, bwaDir, bwaIndex, nBWA, oBWA, tBWA, mBWA, picardDir, saveFiles)
    chrList <- dnaseObject[['chrList']]
    chrom.ref <- dnaseObject[['chrom.ref']]
    dnaseFormat <- tolower(sub(".*\\.", "", dnaseFile))
    if( dnaseFormat %in% c("sam", "bam", "bed") || is.null(bowtieIndex) ){
      dnaseAlign <- list()  ## If DNase file is sam, bam or bed, or align by BWA: No aligment information
    }else{
      dnaseAlign <- .summary(paste(outfileLoc, '/', dnaseName, sep = ""), "priorProcessDNaseBowtie_temp.txt") ##fastq file aligned by Bowtie will have summary information
    }
    
    
    # Process each Histone data     
    link <- vector('list', length(histoneFile))
    names(link) <- histoneName
    histoneAlign <- list()
    for(i in 1:length(histoneFile)){
      link[[i]] <- .histoneProcess(histoneFile[i], fragL, AllocThres, chrList, chrom.ref, capping, outfileLoc = paste(outfileLoc, '/', histoneName[i], '/', sep=''), paste(histoneName[i], '_', outfile, sep=''), bowtieIndex, csemDir, bowtieDir, vBowtie, mBowtie, pBowtie, bwaDir, bwaIndex, nBWA, oBWA, tBWA, mBWA, picardDir, saveFiles)
      
      histoneFormat <- tolower(sub(".*\\.", "", histoneFile[i]))
      if(histoneFormat %in% c("sam", "bam", "bed") || is.null(bowtieIndex)){
        histoneAlign[[histoneName[i]]] <- list()   ## If Histone file is sam, bam or bed, or align using BWA: No aligment information
      }else{
        histoneAlign[[histoneName[i]]] <- .summary(paste(outfileLoc, "/", histoneName[i], sep = ""), "priorProcessHistoneBowtie_temp.txt") #Bowtie alignment has summary information
      }
    }
    
    #generate multipatterns
    s <- scan(chrom.ref, what='char', sep='\n')
    chr_all <- strsplit(s[3], ' ')[[1]]
    chr_length <- as.numeric(strsplit(s[2], ' ')[[1]])
    message( "Info: Partitioning genome based on multiple data sets..." )
    for(i in 1:length(chrList)){
      chr <- chrList[i]
      file_input <- dnaseObject[['posLoc_bychr']][[chr]]    
      for(h in 1:length(histoneFile)){
        file_input <- paste(file_input, ' ', link[[h]][[chr]], sep='')
      }
      # partition genome according to dnase (dnase and histone data)
      script <- 'multi_data_pattern.pl'
      Fn.Path <- system.file(file.path("Perl", script), package="permseq")
      
      CMD<-paste('perl ', Fn.Path, ' ', outfileLoc, '/', chr, '_', outfile, '_positions_cluster.txt', ' ', chr_length[which(chr_all==chr)], ' ', file_input,sep='')
      system(CMD)
    }
    link.t <- vector('list', length(chrList))
    names(link.t) <- chrList
    for(i in chrList){
      link.t[[i]] <- paste(outfileLoc, "/", i, "_", outfile, "_positions_cluster.txt", sep="")
    }
    
    return(new("Prior",
               dnaseName = dnaseName,                 
               dnaseKnots = dnaseObject[['dnaseKnots']],             
               dnaseThres = dnaseObject[['dnaseThres']],                 
               posLoc_bychr = link.t,             
               dataNum = 1 + length(histoneFile),                 
               dnaseAlign = dnaseAlign,
               chrList = chrList,
               histoneAlign = histoneAlign,                
               histoneNum = histoneNum,                 
               histoneName = histoneName,                
               fragL = fragL,                 
               bowtieInfo = bowtieInfo,
               bwaInfo = bwaInfo,
               csemDir = csemDir,
               picardDir = picardDir,
               outfileLoc = outfileLoc,             
               chrom.ref = chrom.ref))
      
  }else if(!is.null(histoneFile) & is.null(dnaseFile)){
    print("use priorHistone_init and priorHistone_multi functions to process Histone data")
  }else{
    print("Provide Files for deriving Priors or run readAllocate directly without prior information.")
  }
  
}
