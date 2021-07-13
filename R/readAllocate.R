readAllocate = function(object = NULL, outfileLoc = "./", outputFormat = NULL, outputFilter = TRUE, chipThres = NULL, chipFile = NULL, bowtieDir = NULL, bowtieIndex = NULL, vBowtie = 2, mBowtie = 99, pBowtie = 8, bwaDir = NULL, bwaIndex = NULL, nBWA = 2, oBWA = 1, tBWA = 8, mBWA = 99, csemDir = NULL, picardDir = NULL, fragL = 200){

  #Refine the input parameters -- NULL NA numeric ""
  #paraList <- c(object, outfileLoc, outputFormat, chipThres, chipFile, bowtieDir, bowtieIndex, vBowtie, mBowtie, pBowtie, bwaDir, bwaIndex, nBWA, oBWA, tBWA, mBWA, csemDir, fragL)
  #for(para in paraList){
  #  para <- .refineParameters(para)
  #}

  #Create the outfileLoc if not exists
    if(!dir.exists(outfileLoc))
        system(paste('mkdir ', outfileLoc, sep=''))
  
    setwd(outfileLoc)
  
  
  # Check parameters
    if(prod(tolower(sub(".*\\.", "", chipFile)) %in% c("fastq", "sam", "bam")) == 0)
        stop("All ChIP-seq files must be in either fastq format or sam , bam format.")

    if(!(tolower(outputFormat) %in% c("tagalign", "bed")))
        stop("outputFormat: Not supported format. readAllocate will generate bam output file, and transform into tagAlign or BED format if specified by outputFormat.")

    
    if(!is.null(chipThres))
        if(chipThres < 0 | chipThres > 1000)
            stop("chipThres must be a number between 0 and 1000.")
  
    if(!is.null(fragL))
        if(fragL<0)
            stop("fragL (fragment length) must be a positive value.")
  
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
    
    
    chipName <- NULL
    if(is.null(object)){ # If object=NULL: Analyze ChIP data without using prior files

        # Warnings
        if(is.null(chipFile))
            stop("Please provide proper path to ChIP-seq file.")
        
        if(sum(tolower(sub(".*\\.", "", chipFile)) == "fastq") > 0 & is.null(bowtieIndex) & is.null(bwaIndex))
            stop("Please provide either Bowtie or BWA index for alignment.")

        if(sum(tolower(sub(".*\\.", "", chipFile)) == "fastq") > 0 & is.null(csemDir))
            stop("Please provide CSEM path for multi-mapping reads allocation.")

        if(sum(tolower(sub(".*\\.", "", chipFile)) == "sam") > 0 & is.null(csemDir))
            stop("Please provide CSEM path for multi-mapping reads allocation.")
        

        for(i in 1:length(chipFile)){#Extract ChIP-seq files names
            t1 <- strsplit(chipFile[i], '/')[[1]]
            t1 <- t1[length(t1)]
            t <- strsplit(t1, paste('.', sub(".*\\.", "", t1), sep=''))[[1]][1]
            chipName <- c(chipName, t)
        }
        chipAlign <- list()
        chipSAM <- c()
        chipBAM <- c()
        for(i in 1:length(chipFile)){
            
            chipFormat <- tolower(sub(".*\\.", "", chipFile[i]))
            # if in bam format: filter by score by request
            if(chipFormat == "bam"){
                print(paste("ChIP-seq file ", chipFile[i], " is in BAM format. Please be sure that this aligned BAM file contains multi-mapping reads already processed. Otherwise, please provide FASTQ format ChIP-seq file.", sep = ""))
                chipAlign[[chipName[i]]] <- NULL
                chipSAM[i] <- NULL
                chipBAM[i] <- chipFile[i]
            }else{
                # if in fastq format: align - remove duplicates - csem - bam

                if(chipFormat == "fastq"){
                    if(!is.null(bowtieIndex)){
                        message( "Info: Aligning reads using Bowtie..." )
                        system(paste(bowtieDir, '/bowtie -q -v ', vBowtie, ' -a -m ', mBowtie, ' -p ', pBowtie, ' --sam ', bowtieIndex, ' ', chipFile[i], ' ', chipName[i], '.sam', " 2>&1 | tee ", outfileLoc, "/priorProcessChIP_", chipName[i], "_Bowtie_temp.txt", sep=''))
                        chipAlign[[chipName[i]]] <- .summary(outfileLoc, paste("priorProcessChIP_", chipName[i], "_Bowtie_temp.txt", sep = ""))
                        chipSAM[i] <- paste(outfileLoc, '/', chipName[i], '.sam', sep='')

                    }else{
                        message( "Info: Aligning reads using BWA..." )
                        system(paste(bwaDir, '/bwa aln -n ', nBWA, ' -o ', oBWA, ' -t ', tBWA, ' ', bwaIndex, ' ', chipFile[i], ' >', chipName[i], '.sai', sep=''))
                        system(paste(bwaDir, '/bwa samse -n ', mBWA, ' ', bwaIndex, ' ', chipName[i], '.sai', ' ', chipFile[i], '  >', chipName[i], '.sam.multiOneLine', sep=''))

                        system(cat("awk '{split($0, tag,\"XA\");split($1, head, \"\"); if (head[1] == \"@\") print $0; else if ($5 >24) print $0; else if (tag[2] != \"\") print $0;}' ", chipName[i], ".sam.multiOneLine | ", bwaDir, "/xa2multi.pl >", chipName[i], ".sam.tmp\n", sep = ""))
                     
                        system(cat("cat ", chipName[i], ".sam.tmp | awk 'BEGIN {FS=\"\\t\" ; OFS=\"\\t\"} ! /^@/ && $6!=\"*\" { cigar=$6; gsub(\"[0-9]+D\",\"\",cigar); n = split(cigar,vals,\"[A-Z]\"); s = 0; for (i=1;i<=n;i++) s=s+vals[i]; seqlen=length($10) ; if (s!=seqlen) {print $0} }' >", chipName[i], ".sam.tmp.badCIGAR\n", sep = ""))

                        system(cat("if [ $(cat ", chipName[i], ".sam.tmp.badCIGAR | wc -l) -gt 0 ]; then\ncat ", chipName[i], ".sam.tmp  | grep -v -F -f ", chipName[i], ".sam.tmp.badCIGAR >", chipName[i], ".sam \nelse\nmv ", chipName[i], ".sam.tmp ", chipName[i], ".sam\nfi\n", sep = "")) 


                        ## system(paste(bwaDir, '/bwa samse -n ', mBWA, ' ', bwaIndex, ' ', chipName[i], '.sai', ' ', chipFile[i], ' | ', bwaDir, '/xa2multi.pl >', chipName[i], '.sam', sep=''))
                        chipAlign[[chipName[i]]] <- NULL
                        chipSAM[i] <- paste(outfileLoc,'/', chipName[i],'.sam',sep='')
                    }
                }else if(chipFormat == "sam"){
                
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
 
                message( "Info: Allocating multi-reads using CSEM..." )
                if(!is.null(bowtieIndex)){
                    CMD <- paste(csemDir,'/run-csem --sam  -p ', pBowtie, ' ', chipSAM[i], ' ', fragL, ' ', chipName[i], '_permseq', sep='') 
                    system(CMD)
                
                }else{
                
                    CMD <- paste(csemDir,'/run-csem --sam  -p ', tBWA, ' ', chipSAM[i], ' ', fragL, ' ', chipName[i], '_permseq', sep='') 
                    system(CMD)
                }
                
                chipBAM[i] <- paste(outfileLoc,'/', chipName[i],'_permseq.bam',sep='')

            }
        }            
    }else{  # If object!=NULL: Analize ChIP data using prior information


        chipName <- object[['chipName']]
        chipSAM <- object[['chipSAM']]
        bowtieIndex <- object[['bowtieInfo']]$bowtieIndex
        pBowtie <- object[['bowtieInfo']]$pBowtie
        tBWA <- object[['bwaInfo']]$tBWA
        csemDir <- object[['csemDir']]
        fragL <- object[['fragL']]
        chipBAM <- c()
        message( "Info: Allocating multi-reads..." )
        for(i in 1:length(chipName)){
            if(!is.null(bowtieIndex)){
                CMD <- paste(csemDir, '/run-csem --sam  -p ', pBowtie, ' --prior ', object[['prior']][i], ' ', chipSAM[i], ' ', fragL, ' ', chipName[i], '_permseq', sep='')
                system(CMD)
            }else{
                CMD <- paste(csemDir, '/run-csem --sam  -p ', tBWA, ' --prior ', object[['prior']][i], ' ', chipSAM[i], ' ', fragL, ' ', chipName[i], '_permseq', sep='')
                system(CMD)

            }
            chipBAM[i] <- paste(outfileLoc,'/',chipName[i],'_permseq.bam',sep='')
        }
        
    }

    # Filter aligned bam file by posterior probability
    if(isTRUE(outputFilter)){
        for(i in 1:length(chipName)){
            bamName <- gsub("\\.bam", "", chipBAM[i])
            # filter low probablity alignments
            if(!is.null(chipThres)){

                CMD <- paste(csemDir, "/sam/samtools view -H ", chipBAM[i], " >tmp.aligned.header", sep = "")
                system(CMD)
                CMD <- paste(csemDir, "/sam/samtools view -F 4 ", chipBAM[i], " | awk -v thres=", chipThres," '{split($NF, zw, \":\");if(zw[3]*1000 > thres) {print $0}}' >tmp.aligned", sep = "")
                print(CMD)
                system(CMD)
                CMD <- "cat tmp.aligned.header tmp.aligned >tmp.alignedFilter.sam"
                system(CMD)
                CMD <- paste(csemDir, "/sam/samtools view -Sb tmp.alignedFilter.sam >", bamName, ".filter.bam", sep="")
                system(CMD)
                CMD <- "rm -rf *tmp.aligned*"
                system(CMD)
            }else{
            # filter to get aligned reads
                CMD <- paste(csemDir, "/sam/samtools view -F 4 -b ", chipBAM[i], " >", bamName, ".filter.bam", sep = "")
                system(CMD)
            }
            
        }
    }

    # Transform aligned bam file to other format
    chipAlignFormat <- NULL
    if(!is.null(outputFormat)){
        message( "Info: Converting BAM to other formats..." )
        if(tolower(outputFormat)=='tagalign'){

            for(i in 1:length(chipName)){
                bamName <- gsub("\\.bam", "", chipBAM[i])
                if(isTRUE(outputFilter)){
                    CMD <- paste(csemDir, '/csem-generate-input --tag-align ', bamName, '.filter.bam', ' ', bamName, '.filter', sep='')
                    system(CMD)
                    chipAlignFormat='.filter.tagAlign'                
                }else{
                    CMD <- paste(csemDir, '/csem-generate-input --tag-align ', bamName, '.bam', ' ', bamName, sep='')
                    system(CMD)
                    chipAlignFormat='.tagAlign'
                }
            }
        }else if(tolower(outputFormat)=='bed'){

            for(i in 1:length(chipName)){
                bamName <- gsub("\\.bam", "", chipBAM[i])
                if(isTRUE(outputFilter)){
                    CMD <- paste(csemDir, '/csem-generate-input --bed ', bamName, '.filter.bam', ' ', bamName, '.filter', sep='')
                    system(CMD)
                    chipAlignFormat='.filter.bed'
                }else{
                    CMD <- paste(csemDir, '/csem-generate-input --bed ', bamName, '.bam', ' ', bamName, sep='')
                    system(CMD)
                    chipAlignFormat='.bed'
                }
            }
        }else{
            print('Not supported format! Please choose between tagAlign and BED.')
        }
        chipAlignFormat <- paste(gsub("\\.bam", "", chipBAM), chipAlignFormat, sep='')
    }
    
  # Summarize information
    if(length(object) == 0){
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
      
        return(new(Class = "Prior",
                   chipAlign = chipAlign,
                   chipSAM = chipSAM,
                   chipAllocate = paste(outfileLoc, "/", chipName, '_permseq.bam', sep=''),
                   chipAlignFormat = paste(outfileLoc, "/", chipAlignFormat, sep = ""),
                   outfileLoc = outfileLoc,
                   chipName = chipName,
                   chipNum = length(chipFile),
                   fragL = fragL,
                   bowtieInfo = bowtieInfo,
                   bwaInfo = bwaInfo,
                   csemDir = csemDir,
                   picardDir = picardDir,
                   dataNum = length(chipFile)))
        
    }else{

        dnaseThres <- object@dnaseThres
        dnaseKnots <- object@dnaseKnots
        
        outfile_chipmean <- paste("chip", 1:length(chipName), "_chipmean", sep="")
        
        posLoc_bychr <- vector('list', length(object@dnaseName))
        names(posLoc_bychr) <- object@dnaseName
        posLoc_bychr[[length(object@dnaseName)]] <- object@posLoc_bychr
        
    #calculate averaged chip read counts according to object (dnase or histone information) so that we can save time for plotting.
        if(!file.exists((paste(outfileLoc, '/', names(posLoc_bychr)[1], '_', outfile_chipmean[1], sep='')))){
            .chipMeanCounts(object, posLoc_bychr, chipSAM, paste(chipName, '.sam', sep=''), outfileLoc, outfile_chipmean)
        }
        
        object['chipAllocate'] <- paste(outfileLoc, "/", chipName, '_permseq.bam', sep='')
        object['chipAlignFormat'] <- paste(outfileLoc, "/", chipAlignFormat, sep = "")
        object['outfileLoc'] <- outfileLoc
        
        return(object)
    }
    message( "Info: ---- Done! ----" )
}

