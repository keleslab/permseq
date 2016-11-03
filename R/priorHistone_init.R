priorHistone_init = function(histoneFile = NULL, histoneName = NULL, chipFile = NULL, fragL = 200, AllocThres = 900, chrList = NULL, capping = 0, outfileLoc = "./",  bowtieDir = NULL, bowtieIndex = NULL,  vBowtie = 2, mBowtie = 99, pBowtie = 8, bwaDir = NULL, bwaIndex = NULL, nBWA = 2, oBWA = 1, tBWA = 8, mBWA = 99, csemDir = NULL, saveFiles = TRUE){

  if(length(chipFile)>1){
    message( "More than one ChIP-seq data detected: Replicates must be from the same TF (experiment). Otherwise, we recommend run different ChIP-seq data separtely.")
  }else if(is.null(chipFile)){
    stop("Please provide proper path to ChIP-seq file.")
  }
  
  # Check parameters
  if(fragL <= 0)
    stop("fragL (fragment length) must be a positive value.")

  if(AllocThres<0 | AllocThres>1000)
    stop("AllocThres (allocation score threshold) must be a number between 0 and 1000.")

  if(is.null(histoneFile))
    stop("Please provide proper path to histone file or use priorProcess for prior building using only DNase-seq or readAllocate for no prior Multi-reads allocation.")

  if(prod(tolower(sub(".*\\.", "", histoneFile)) %in% c("fastq", "sam", "bam", "bed")) == 0 )
    stop("All Histone ChIP-seq files must be in fastq, sam , bam or bed format.")

  if(sum(tolower(sub(".*\\.", "", histoneFile)) == "fastq") > 0 & is.null(bowtieIndex) & is.null(bwaIndex))
      stop("Please provide either Bowtie or BWA index for alignment.")
  
  if(prod(tolower(sub(".*\\.", "", chipFile)) %in% c("fastq", "sam")) == 0)
    stop("All ChIP-seq files must be in either fastq format or sam format.")

  if(sum(tolower(sub(".*\\.", "", chipFile)) == "fastq") > 0 & is.null(bowtieIndex) & is.null(bwaIndex))
      stop("Please provide either Bowtie or BWA index for alignment.")



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
  
  #If the output folder does not exist, we create it automatically
  if(!file.exists(outfileLoc)){
    system(paste('mkdir ', outfileLoc, sep=''))
  }
  #change into the output folder directory
  setwd(outfileLoc)
  
  #build the chromosome reference file using bowtie-inspect
  if(!is.null(bowtieIndex)){
    script <- "genRef.pl"
    Fn.Path <- system.file(file.path("Perl", script), package = "permseq")
    bowtieIndexName <-  gsub(".*/(.*)", "\\1", bowtieIndex)
    chrom.ref <- paste(outfileLoc, "/", bowtieIndexName, ".ref", sep = "")
    CMD <- paste(bowtieDir, "/bowtie-inspect -s ", bowtieIndex, " | perl ", Fn.Path, " ", chrom.ref, sep = "")
    system(CMD, intern = TRUE)
  }else{
    chrom.ref <- NULL
  }
  
  # Sumarize Bowtie/BWA Input Parameters Information
  if(is.null(bwaIndex)){
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
  
  
  #If the histone name(s) is not given, use the histone directory name
  if(is.null(histoneName)){
    histoneName <- 1:length(histoneFile)
  }
  link <- vector('list', length(histoneFile))
  names(link) <- histoneName
  histoneAlign <- list()

  # Process each Histone data - align to the reference genome and give the bowtie summary information
  for(i in 1:length(histoneFile)){
    link[[i]] <- .multihistoneProcess(histoneFile[i], fragL, AllocThres, chrList, chrom.ref, capping=0, outfileLoc = paste(outfileLoc, '/', histoneName[i], '/', sep=''), histoneName[i], bowtieIndex, csemDir, bowtieDir, vBowtie, mBowtie, pBowtie, bwaDir, bwaIndex, nBWA, oBWA, tBWA, mBWA, saveFiles)
    chrList <- link[[i]][['chrList']]
    histoneFormat <- tolower(sub(".*\\.", "", histoneFile[i]))
    if(is.null(bowtieIndex) || histoneFormat %in% c("sam", "bam", "bed")){ 
      histoneAlign[[histoneName[i]]] <- list()  ## If Histone file is bam or bed: No aligment information
    }else{
      histoneAlign[[histoneName[i]]] <- .summary(paste(outfileLoc, "/", histoneName[i], sep = ""), "priorProcessHistoneBowtie_temp.txt")
    }
    
  }
  
  # Align each ChIP data
  outfile_chip <- NULL
  chipAlign <- list()
  chipSAM <- c()
  
  #Get the ChIP-seq data name
  for(i in 1:length(chipFile)){
    t1 <- strsplit(chipFile[i], '/')[[1]]
    t1 <- t1[length(t1)]
    t <- strsplit(t1, paste('.', sub(".*\\.", "", t1), sep=''))[[1]][1]
    outfile_chip <- c(outfile_chip, t)
  }

  if(!is.null(bowtieIndex)){
    message("Info: Align reads using Bowtie...")
    for(i in 1:length(chipFile)){
      chipFormat <- tolower(sub(".*\\.", "", chipFile[i]))
      if(chipFormat == "fastq"){
        # Bowtie alignment
        system(paste(bowtieDir, '/bowtie -q -v ', vBowtie, ' -a -m ', mBowtie, ' -p ', pBowtie, ' --sam ', bowtieIndex, ' ', chipFile[i], ' ', outfile_chip[i], '.sam',  " 2>&1 | tee ", outfileLoc, "/priorProcessChIP_", outfile_chip[i], "_Bowtie_temp.txt", sep=''))
    
        # Summarize ChIP alignment information 
        chipAlign[[outfile_chip[i]]] <- .summary(outfileLoc, paste("priorProcessChIP_", outfile_chip[i], "_Bowtie_temp.txt", sep = ""))
        chipSAM[i] <- paste(outfileLoc, "/", outfile_chip[i], ".sam", sep = "")
      }else{
        print(paste("ChIP-seq file ", outfile_chip[i], " is in SAM format. The preprocessed alignment should contain multi-mapping reads. Otherwise, please provide FASTQ format ChIP-seq file.", sep = ""))
        
        chipAlign[[outfile_chip[i]]] <- NULL
        chipSAM[i] <- chipFile[i]
      }
    }
  }else{
    message( "Info: Aligning reads using BWA..." )
      
    for(i in 1:length(chipFile)){
      chipFormat <- tolower(sub(".*\\.", "", chipFile[i]))
      if(chipFormat == "fastq"){
      
        system(paste(bwaDir, '/bwa aln -n ', nBWA, ' -o ', oBWA, ' -t ', tBWA, ' ', bwaIndex, ' ', chipFile[i], ' >', outfile_chip[i], '.sai', sep=''))
        system(paste(bwaDir, '/bwa samse -n ', mBWA, ' ', bwaIndex, ' ', outfile_chip[i], '.sai', ' ', chipFile[i], ' | ', bwaDir, '/xa2multi.pl >', outfile_chip[i], '.sam', sep=''))
        
        chipSAM[i] <- paste(outfileLoc,'/',outfile_chip[i],'.sam',sep='')
         }else{
           
           print(paste("ChIP-seq file ", outfile_chip[i], " is in SAM format. The preprocessed aligned SAM file should contain multi-mapping reads. Otherwise, please provide FASTQ format ChIP-seq file.", sep = ""))
           
           chipSAM[i] <- chipFile[i]
         }
         chipAlign[[outfile_chip[i]]] <- NULL
    }
  }
  
  outfile_chipmean <- paste(outfile_chip, '_chipmean_temp', sep='')

  # Calcualte averaged chip counts at different groups
  link.temp <- vector('list', length(histoneName))
  names(link.temp) <- histoneName
  for(i in 1:length(histoneName)){
    link.temp[[i]] <- link[[i]][['posLoc_bychr']]
  }
  .chipMeanCounts(prior = NULL, link.temp, chipSAM, paste(outfile_chip, '.sam', sep=''), outfileLoc, outfile_chipmean)

  #prepare for the plot to select the best histone to be used as dnase data
  ylim <- vector('list', length(chipFile))
  for(i in 1:length(chipFile))
    ylim[[i]] <- c(0, 0)
  reps <- vector('list', length(histoneFile))
  for(i in 1:length(reps)){
    reps[[i]] <- vector('list', length(chipFile))
    for(j in 1:length(chipFile)){
      reps[[i]][[j]] <- read.table(paste(histoneName[i], '_', outfile_chipmean[j], sep=''))
      ylim[[j]][1] <- min(c(ylim[[j]][1], min(unlist(reps[[i]][[j]][2, ]/reps[[i]][[j]][1, ]))))
      ylim[[j]][2] <- max(c(ylim[[j]][2], max(unlist(reps[[i]][[j]][2, ]/reps[[i]][[j]][1, ]))))  
    }
  }
  # Generate plot used to select one histone data which will be used as "DNase data" in the following step.
  pdf(paste(outfileLoc, '/', 'marginal_histone_plot.pdf', sep=''))
  for(j in 1:length(outfile_chip))
    for(i in 1:length(reps)){
      .fitPlot(list(reps[[i]][[j]]), link[[i]][['histoneThres']], link[[i]][['histoneKnots']],name = outfile_chip[j], xlab = histoneName[i], ylim=ylim[[j]])
    }
  dev.off()
  
  # Summarizing information and create "Prior" class. The information saved in the "Prior" object will be passed to "priorHistone_multi"
  chipUni <- c()
  for(i in 1:length(chipFile)){
   chipUni[i] <-  paste(outfileLoc, '/', outfile_chip[i], '.sam.uni.bed', sep='')
  }

  #result <- vector('list', length(histoneName))
  #names(result) <- histoneName
  #for(i in 1:length(result)){
  #  result[[i]] <- vector('list', length(outfile_chip))
  #  names(result[[i]]) <- outfile_chip
  #  for(j in 1:length(chipFile)){
  #    result[[i]][[j]] <- unlist(reps[[i]][[j]][2,]/reps[[i]][[j]][1,])
  #  }
  #}
  #link[['meanchip_histone']] <- result
 
  system(paste('rm -rf ', outfileLoc, '/*temp', sep=''))

  #create new "Prior" object
  result = new("Prior",
    chipName = outfile_chip,
    chipSAM = chipSAM,
    chipUni = chipUni,
    chipNum = length(chipFile),
    histoneName = histoneName,
    fragL = fragL,
    chrList = chrList,
    bowtieInfo = bowtieInfo,
    bwaInfo = bwaInfo,
    csemDir = csemDir,
    chrom.ref = chrom.ref,
    histoneAlign = histoneAlign,
    chipAlign = chipAlign,
    histoneNum = length(histoneFile),
    outfileLoc = outfileLoc,
    dataNum = length(histoneName)
    )
  for(i in 1:length(histoneFile)){
    result@dnaseHistone[[histoneName[i]]] <- link[[histoneName[i]]]
  }
  return(result)
}

#plot the histone vs ChIP-seq read counts
.fitPlot = function(reps, histoneThres, histoneKnots, name, xlab, ylim){
  par(mar = c(5, 5, 4, 2) + 0.1)
  for (i in 1:length(reps)) {
    a = reps[[i]][2, ]/reps[[i]][1, ]
    matplot(histoneThres, t(a), col = "black", lwd = 1, pch = 20, ylab = "Average ChIP read count", xlab = paste("Histone (", xlab, ") read count", sep = ""), main=name, ylim=ylim)
    for (j in 1:length(histoneKnots)) abline(v = histoneKnots[j], col = "blue", lwd = 2)
    }
}
