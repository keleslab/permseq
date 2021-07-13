priorGenerate = function(object = NULL, chipFile = NULL, maxHistone = 6, outfileLoc="./"){

  #Refine the input parameters -- NULL NA numeric ""
  #paraList <- c(object, chipFile, maxHistone, outfileLoc)
  #for(para in paraList){
  #  para <- .refineParameters(para)
  #}
  
  #Check parameters
    if(is.null(chipFile))
        stop("Please provide the right path to ChIP-seq file.")
    
    if(length(chipFile) > 1)
        message( "More than one ChIP-seq data detected: Replicates must be from the same TF (experiment). Otherwise, we recommend run different ChIP-seq data separately.")

    if(prod(tolower(sub(".*\\.", "", chipFile)) %in% c("fastq", "sam")) == 0)
        stop("All ChIP-seq files must be in either fastq format or sam format. If it is aligned BAM file with permseq or CSEM aligned multi-mapping reads, readAllocate() function can be employed to transform BAM to other format say tagAlign or BED.")
    
    if(is.null(object))
        stop("Please provide the returned object from priorProcess step.")

  
    csemDir <- object['csemDir']
    picardDir <- object['picardDir']
    chrlist <- names(object['posLoc_bychr'])
    dnaseKnots <- object['dnaseKnots']
    dnaseThres <- object['dnaseThres']
    object['chipNum'] <- length(chipFile)
    bwaInfo <- object['bwaInfo']
  
    #If the output folder does not exist, we create it automatically
    if(!dir.exists(outfileLoc)){
        system(paste('mkdir ', outfileLoc, sep=''))
    }
    #change into the output folder directory
    setwd(outfileLoc)
    
    chipName <- NULL
    # Name of ChIP file(s)
    for(i in 1:length(chipFile)){
        t1 <- strsplit(chipFile[i],'/',)[[1]]
        t1 <- t1[length(t1)]
        t <- strsplit(t1,paste('.',sub(".*\\.","",t1), sep=''))[[1]][1]
        chipName <- c(chipName, t)
    }
  
    if(length(object['chipSAM'])==0){ # If ChIP data hasn't been processed before
        
        chipAlign <- list()
        chipSAM <- c()
        if(is.null(bwaInfo)){
            message( "Info: Aligning reads using Bowtie..." )      
            # Obtain information from "object"
            bowtieIndex <- object['bowtieInfo']$bowtieIndex
            bowtieDir <- object['bowtieInfo']$bowtieDir
            vBowtie <- object['bowtieInfo']$vBowtie
            mBowtie <- object['bowtieInfo']$mBowtie
            pBowtie <- object['bowtieInfo']$pBowtie
            for(i in 1:length(chipFile)){
                chipFormat <- tolower(sub(".*\\.", "", chipFile[i]))
                if(chipFormat == "fastq"){
                    system(paste(bowtieDir, '/bowtie -q -v ', vBowtie, ' -a -m ', mBowtie, ' -p ', pBowtie, ' --sam ', bowtieIndex, ' ', chipFile[i], ' ', chipName[i], '.sam', " 2>&1 | tee ", outfileLoc, "/priorProcessChIP_", chipName[i], "_Bowtie_temp.txt", sep=''))
                    chipAlign[[chipName[i]]] <- .summary(outfileLoc, paste("priorProcessChIP_", chipName[i], "_Bowtie_temp.txt", sep = ""))
                    chipSAM[i] <- paste(outfileLoc,'/', chipName[i],'.sam',sep='')
                }else{
                    
                    print(paste("ChIP-seq file ", chipFile[i], " is in SAM format. The preprocessed alignment should contain multi-mapping reads. Otherwise, please provide FASTQ format ChIP-seq file.", sep = ""))
        
                    chipAlign[[chipName[i]]] <- NULL
                    chipSAM[i] <- chipFile[i]
                }

                
                ## Remove duplicates from SAM file
                file.sam <- chipSAM[i]
                filename <- gsub("\\.sam", "", chipSAM[i])
                if(!is.null(picardDir)){
                    CMD <- "mkdir -p TMP_permseq"
                    system(CMD)
                    CMD <- paste("java -Djava.io.tmpdir=TMP_permseq -jar ", picardDir, " SortSam I=", file.sam, " O=", filename, ".sort.sam VALIDATION_STRINGENCY=LENIENT SORT_ORDER=coordinate TMP_DIR=TMP_permseq", sep = "")
                    system(CMD)
                    CMD <- paste("java -Djava.io.tmpdir=TMP_permseq -jar ", picardDir, " MarkDuplicates I=", filename, ".sort.sam ", "O=", filename, ".nodup.sam METRICS_FILE=QCmetric VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=true TMP_DIR=TMP_permseq", sep = "")
                    system(CMD)
                    CMD <- paste("java -Djava.io.tmpdir=TMP_permseq -jar ", picardDir, " SortSam I=", filename, ".nodup.sam O=", filename, ".nodup.sort.sam VALIDATION_STRINGENCY=LENIENT SORT_ORDER=queryname TMP_DIR=TMP_permseq", sep = "")
                    system(CMD)
                    CMD <- "rm -rf TMP_permseq"
                    system(CMD)

                }else{
                    CMD <- paste(csemDir, "/sam/samtools view -Sb -h ", file.sam ," >", filename, ".tmp.bam", sep = "")
                    system(CMD)
                    CMD <- paste(csemDir, "/sam/samtools sort ", filename, ".tmp.bam ", filename, ".tmp.sort", sep = "")
                    system(CMD)
                    
                    CMD <- paste(csemDir, "/sam/samtools rmdup ", filename, ".tmp.sort.bam ", filename, ".tmp.nodup.bam", sep = "")
                    system(CMD)
                    CMD <- paste(csemDir, "/sam/samtools view -h ", filename, ".tmp.nodup.bam >", filename, ".tmp.nodup.sam", sep = "")
                    system(CMD)
                    CMD <- paste(csemDir, "/sam/samtools sort -n ", filename, ".tmp.nodup.sam ", filename, ".nodup.sort", sep = "")
                    system(CMD)
                
#                CMD <- paste("rm -rf ", filename, ".tmp*",  sep = "")
#                system(CMD)
                }
                chipSAM[i] <- paste(filename, ".nodup.sort.sam", sep = "")
            
                

            }
        }else{
            message( "Info: Aligning reads using BWA..." )
            bwaIndex <- object['bwaInfo']$bwaIndex
            bwaDir <- object['bwaInfo']$bwaDir
            nBWA <- object['bwaInfo']$nBWA
            oBWA <- object['bwaInfo']$oBWA
            tBWA <- object['bwaInfo']$tBWA
            mBWA <- object['bwaInfo']$mBWA
            
            for(i in 1:length(chipFile)){
                chipFormat <- tolower(sub(".*\\.", "", chipFile[i]))
                if(chipFormat == "fastq"){
                    
                    system(paste(bwaDir, '/bwa aln -n ', nBWA, ' -o ', oBWA, ' -t ', tBWA, ' ', bwaIndex, ' ', chipFile[i], ' >', chipName[i], '.sai', sep=''))
                    system(paste(bwaDir, '/bwa samse -n ', mBWA, ' ', bwaIndex, ' ', chipName[i], '.sai', ' ', chipFile[i], '  >', chipName[i], '.sam.multiOneLine', sep=''))

                    system(cat("awk '{split($0, tag,\"XA\");split($1, head, \"\"); if (head[1] == \"@\") print $0; else if ($5 >24) print $0; else if (tag[2] != \"\") print $0;}' ", chipName[i], ".sam.multiOneLine | ", bwaDir, "/xa2multi.pl >", chipName[i], ".sam.tmp\n", sep = ""))
                     
                    system(cat("cat ", chipName[i], ".sam.tmp | awk 'BEGIN {FS=\"\\t\" ; OFS=\"\\t\"} ! /^@/ && $6!=\"*\" { cigar=$6; gsub(\"[0-9]+D\",\"\",cigar); n = split(cigar,vals,\"[A-Z]\"); s = 0; for (i=1;i<=n;i++) s=s+vals[i]; seqlen=length($10) ; if (s!=seqlen) {print $0} }' >", chipName[i], ".sam.tmp.badCIGAR\n", sep = ""))

                    system(cat("if [ $(cat ", chipName[i], ".sam.tmp.badCIGAR | wc -l) -gt 0 ]; then\ncat ", chipName[i], ".sam.tmp  | grep -v -F -f ", chipName[i], ".sam.tmp.badCIGAR >", chipName[i], ".sam \nelse\nmv ", chipName[i], ".sam.tmp ", chipName[i], ".sam\nfi\n", sep = "")) 

 
                    
                    ## system(paste(bwaDir, '/bwa samse -n ', mBWA, ' ', bwaIndex, ' ', chipName[i], '.sai', ' ', chipFile[i], ' | ', bwaDir, '/xa2multi.pl >', chipName[i], '.sam', sep=''))
                    
                    chipSAM[i] <- paste(outfileLoc,'/', chipName[i],'.sam',sep='')
                }else{
                    print(paste("ChIP-seq file ", chipFile[i], " is in SAM format. The preprocessed aligned SAM file should contain multi-mapping reads. Otherwise, please provide FASTQ format ChIP-seq file.", sep = ""))
        
                    chipSAM[i] <- chipFile[i]
                }
                chipAlign[[chipName[i]]] <- NULL

                ## Remove duplicates from SAM file
                file.sam <- chipSAM[i]
                filename <- gsub("\\.sam", "", chipSAM[i])
                if(!is.null(picardDir)){
                    CMD <- "mkdir -p TMP_permseq"
                    system(CMD)
                    CMD <- paste("java -Djava.io.tmpdir=TMP_permseq -jar ", picardDir, " SortSam I=", file.sam, " O=", filename, ".sort.sam VALIDATION_STRINGENCY=LENIENT SORT_ORDER=coordinate TMP_DIR=TMP_permseq", sep = "")
                    system(CMD)
                    CMD <- paste("java -Djava.io.tmpdir=TMP_permseq -jar ", picardDir, " MarkDuplicates I=", filename, ".sort.sam ", "O=", filename, ".nodup.sam METRICS_FILE=QCmetric VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=true TMP_DIR=TMP_permseq", sep = "")
                    system(CMD)
                    CMD <- paste("java -Djava.io.tmpdir=TMP_permseq -jar ", picardDir, " SortSam I=", filename, ".nodup.sam O=", filename, ".nodup.sort.sam VALIDATION_STRINGENCY=LENIENT SORT_ORDER=queryname TMP_DIR=TMP_permseq", sep = "")
                    system(CMD)
                    CMD <- "rm -rf TMP_permseq"
                    system(CMD)

                }else{
                    CMD <- paste(csemDir, "/sam/samtools view -Sb -h ", file.sam ," >", filename, ".tmp.bam", sep = "")
                    system(CMD)
                    CMD <- paste(csemDir, "/sam/samtools sort ", filename, ".tmp.bam ", filename, ".tmp.sort", sep = "")
                    system(CMD)
                    
                    CMD <- paste(csemDir, "/sam/samtools rmdup ", filename, ".tmp.sort.bam ", filename, ".tmp.nodup.bam", sep = "")
                    system(CMD)
                    CMD <- paste(csemDir, "/sam/samtools view -h ", filename, ".tmp.nodup.bam >", filename, ".tmp.nodup.sam", sep = "")
                    system(CMD)
                    CMD <- paste(csemDir, "/sam/samtools sort -n ", filename, ".tmp.nodup.sam ", filename, ".nodup.sort", sep = "")
                    system(CMD)
                
#                CMD <- paste("rm -rf ", filename, ".tmp*",  sep = "")
#                system(CMD)
                }
                chipSAM[i] <- paste(filename, ".nodup.sort.sam", sep = "")
           
            }
            
        }
        object['chipAlign'] <- chipAlign
        object['chipSAM'] <- chipSAM
    }


    outfile_chipmean <- paste(chipName, '_chipmean_temp',sep='')
  
    message( "Info: Calculating averaged ChIP counts..." )
    object <- .chipMeanCounts_multi(object, outfileLoc, outfile_chipmean)
  
    message( "Info: Variable selection..." )
    num_require <- 2 + length(dnaseKnots) + maxHistone*2
    if(object['dataNum'] > 1){ 
        setwd(outfileLoc)
        index <- c(rep(1, 4), NA, rep(2:object['dataNum'], each=2))
         # Group Lasso
        result <- mclapply(outfile_chipmean, function(x){.run_grplasso(x, dnaseKnots, dnaseThres, chrlist, index, num_require)}, mc.preschedule = FALSE)
        choose <- data_all <- vector('list', length(outfile_chipmean))
        
        for(i in 1:length(outfile_chipmean)){
         # For each ChIP data - choose: selected histone(s) data, data_all_combined: aggregated data
            choose[[i]] <- result[[i]][['choose']]
            data_all[[i]] <- result[[i]][['data_all_combined']]
        }
        
        inter <- filter <- vector('list', min(unlist(lapply(choose, length))))
        k <- min(c(object['dataNum']-1, maxHistone))
        inter <- choose[[1]][1:k]
        if(length(object['chipSAM']) > 1){
            for(j in 2:(length(object['chipSAM'])))
                inter <- unique(intersect(inter, choose[[j]][1:k]))
        }
        inter <- inter[order(as.numeric(inter), decreasing=F)]
        if(length(inter) != 0){
            filter <- c(1:(object['dataNum']-1))[-inter]
        }else{
            filter <- c(1:(object['dataNum']-1))
            inter <- NULL
        }
        if(length(filter) == 0)
            filter=NULL
    }else{
        data_all <- vector('list', length(outfile_chipmean))
        for(i in 1:length(outfile_chipmean)){
            file <- outfile_chipmean[i]
            for(j in 1:length(chrlist)){
                chr.i <- chrlist[j]
                data_all[[i]] <- rbind(data_all[[i]], t(.fast.read.table(paste(outfileLoc, '/', chr.i,'_', file, sep=''), 3)))
            }
            data_all[[i]] <- .aggregate_data(data_all[[i]])
        }
        filter <- NULL
        inter <- NULL
    }
    
    message( "Info: Model fitting..." )
    y <- vector('list', length(chipFile))
    
    for(i in 1:length(chipFile)){
        y[[i]] <- .fitRlm(data_all[[i]], filter, object)
    }

    message( "Info: Generating prior files..." )
    setwd(outfileLoc)
    CMD <- 'cat '
    for(i in 1:length(object['posLoc_bychr'])){
        CMD <- paste(CMD, object['posLoc_bychr'][[i]], sep=' ')
    }
    CMD <- paste(CMD, '>cluster_temp.txt', sep='')
    system(CMD)
    
    positionFileLoc_tfs <- 'cluster_temp.txt'
    tt <- 0:(length(object['dnaseThres'])-1)
    tt <- tt[order(as.numeric(tt), decreasing=F)]
    if(object['dataNum']>1){ 
        data_id <- t(apply(as.matrix(data_all[[1]][, 1]), 1, function(x){strsplit(as.character(x), '_')[[1]]}))
        tt <- data_id[, 1]
        for(i in inter){
            tt <- paste(tt, data_id[, 1+i], sep='_')
        }
        tt <- unique(tt)
        tt <- matrix(as.character(tt), nrow=1)
        write.table(tt, file='dnase_pattern_id_temp.txt', quote=F, col.names=F, row.names=F, sep='\t')
        script <- 'select_multi_pattern.pl'
        Fn.Path <- system.file( file.path("Perl", script), package="permseq")
        CMD <- paste('perl  ', Fn.Path, 'dnase_pattern_id_temp.txt', 'cluster_temp2.txt', positionFileLoc_tfs, sep=' ')
        for(i in c(0,inter)){
            CMD <- paste(CMD, i, sep=' ')
        }
        system(CMD)
        positionFileLoc_tfs <- 'cluster_temp2.txt'
    }
    for(i in 1:length(chipFile)){
        parameter_a <- c(length(tt), as.numeric(y[[i]][match(tt, y[[i]][, 1]), 2])+1) 
        write.table(matrix(parameter_a,nrow=1), file='rep_dnase_parameter_temp.txt', col.names=F, row.names=F, quote=F)
        system(paste('cat rep_dnase_parameter_temp.txt ', positionFileLoc_tfs, '>prior_', chipName[i], '.txt', sep=''))
    }
    system('rm -rf *temp*')
    

    dnaseThres <- object@dnaseThres
    dnaseKnots <- object@dnaseKnots
  
    outfile_chipmean <- paste("chip", 1:length(chipName), "_chipmean", sep="")
    
    posLoc_bychr <- vector('list', length(object@dnaseName))
    names(posLoc_bychr) <- object@dnaseName
    posLoc_bychr[[length(object@dnaseName)]] <- object@posLoc_bychr
    #calculate averaged chip read counts according to object (dnase or histone information)
    if(!file.exists((paste(outfileLoc, '/', names(posLoc_bychr)[1], '_', outfile_chipmean[1], sep='')))){
        .chipMeanCounts(object, posLoc_bychr, chipSAM, paste(chipName, ".nodup.sam", sep = ""), outfileLoc, outfile_chipmean)
    }
    if(is.null(inter)){
        object['prior'] <- paste(outfileLoc, '/', 'prior_', chipName, '.txt', sep='')
        object['outfileLoc'] <- outfileLoc
        object['chipName'] <- chipName
        
    }else{
        object['prior'] <- paste(outfileLoc, '/', 'prior_', chipName, '.txt', sep='')
        object['outfileLoc'] <- outfileLoc
        object['histoneGrpL'] <- object['histoneName'][inter]
        object['chipName'] <- chipName
    }
    return(object)
}
