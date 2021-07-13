.dnaseProcess = function(dnaseFile, fragL, AllocThres, chrList, chrom.ref, capping, outfileLoc, outfile, bowtieIndex, csemDir, bowtieDir, vBowtie, mBowtie, pBowtie, bwaDir, bwaIndex, nBWA, oBWA, tBWA, mBWA, picardDir, saveFiles){

    if(!file.exists(outfileLoc)){
        system(paste('mkdir ', outfileLoc, sep=''))
    }
    setwd(outfileLoc)
  
    dnaseFormat <- tolower(sub(".*\\.", "", dnaseFile))
  
  # Define filenames and extensions  
    filename <- strsplit(dnaseFile, "/")[[1]]
    filename <- filename[length(filename)]
    filename <- strsplit(filename, paste('.', sub(".*\\.", "", filename), sep=''))[[1]][1]
    file.sai <- paste(filename, '.sai', sep='')
    file.sam <- paste(filename, '.sam', sep='')
    file.bam <- paste(filename, '.bam', sep='')
    file.bed <- paste(filename, '.bed', sep='')
    file.bed_AllocThres <- paste(filename, '_', AllocThres, '.bed', sep='')
  
  # Deal with different DNase-seq file formats (fastq, sam, bam or bed)
    if(dnaseFormat != "bed"){
        if(dnaseFormat != "bam"){
            if(dnaseFormat != "sam"){ #DNase-seq file is in fastq format - align - CSEM allocate multi-reads
                if(!is.null(bowtieIndex)){
                    message( "Info: Aligning DNase-seq reads using Bowtie." )
                    ## Generates the SAM file
                    system(paste(bowtieDir, '/bowtie -q -v ', vBowtie, ' -a -m ', mBowtie, ' -p ', pBowtie, ' --sam ', bowtieIndex, ' ', dnaseFile, ' ', file.sam, " 2>&1 | tee ", outfileLoc, "/priorProcessDNaseBowtie_temp.txt", sep=''))  
                }else{
                    message( "Info: Aligning DNase-seq reads using BWA." )
                    ## Generates the SAM file
                    system(paste(bwaDir, '/bwa aln -n ', nBWA, ' -o ', oBWA, ' -t ', tBWA, ' ', bwaIndex, ' ', dnaseFile, ' >', file.sai, sep=''))
                    system(paste(bwaDir, '/bwa samse -n ', mBWA, ' ', bwaIndex, ' ', outfileLoc, "/", file.sai, ' ', dnaseFile, ' >', file.sam, '.multiOneLine', sep=''))

                    system(cat("awk '{split($0, tag,\"XA\");split($1, head, \"\"); if (head[1] == \"@\") print $0; else if ($5 >24) print $0; else if (tag[2] != \"\") print $0;}' ", file.sam, ".multiOneLine | ", bwaDir, "/xa2multi.pl >", file.sam, ".tmp\n", sep = ""))
                     
                    system(cat("cat ", file.sam, ".tmp | awk 'BEGIN {FS=\"\\t\" ; OFS=\"\\t\"} ! /^@/ && $6!=\"*\" { cigar=$6; gsub(\"[0-9]+D\",\"\",cigar); n = split(cigar,vals,\"[A-Z]\"); s = 0; for (i=1;i<=n;i++) s=s+vals[i]; seqlen=length($10) ; if (s!=seqlen) {print $0} }' >", file.sam, ".tmp.badCIGAR\n", sep = ""))

                    system(cat("if [ $(cat ", file.sam, ".tmp.badCIGAR | wc -l) -gt 0 ]; then\ncat ", file.sam, ".tmp  | grep -v -F -f ", file.sam, ".tmp.badCIGAR >", file.sam, "\nelse\nmv ", file.sam, ".tmp ", file.sam, "\nfi\n", sep = "")) 

                    
                    
                }
            }else{#DNase-seq file is in SAM format - can skip alignment step and directly allocate multi-reads using CSEM
                print(paste("DNase-seq file ", dnaseFile, " is in SAM format. The preprocessed aligned SAM file should contain multi-mapping reads. Otherwise, please provide FASTQ format DNase-seq file.", sep = ""))
                file.sam <- dnaseFile
            }

            ## Get chrList from SAM file
            if(is.null(chrList)){ 
                con <- file(file.sam, open = "r")
                chrList <- c()
                n <- 1
                while(substring(oneLine <- readLines(con, n=1, warn = FALSE), 1, 3) != "@PG"){
                    if(grepl("SN", oneLine)){
                        chrList[n] <- gsub('.*SN:(.*)\t.*', '\\1', oneLine)
                        n <- n+1
                    } 
                }
            }

            ## Remove duplicates from SAM file
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
            file.sam <- paste(filename, ".nodup.sort.sam", sep = "")
            
            message( "Info: Allocating DNase-seq multi-reads using CSEM" )
            
            if(!is.null(bowtieIndex)){
                
                system(paste(csemDir, "/run-csem --sam -p ", pBowtie, ' ', file.sam, ' ', fragL, ' ', filename, sep=''))
                
            }else{

                system(paste(csemDir, "/run-csem --sam -p ", tBWA, ' ', file.sam, ' ', fragL, ' ', filename, sep=''))

            }

        }else{#DNase-seq file is in BAM format - It is better to be the CSEM alignment results.
            print(paste("DNase-seq file ", dnaseFile, " is in BAM format. The preprocessed aligned BAM file should contain multi-mapping reads allocated by CSEM. Otherwise, please provide FASTQ format DNase-seq file.", sep = ""))
            file.bam <- dnaseFile
        }
        message( "Info: Transfering DNase-seq multi-reads alignment bam file into bed file using CSEM" )
        system(paste(csemDir, "/csem-generate-input --bed ", file.bam, ' ', filename, sep=''))
    }else{#DNase-seq file is in BED format - It is better to be the CSEM alignment results.  
        print(paste("DNase-seq file ", dnaseFile, " is in SAM format. The preprocessed aligned BED file should contain multi-mapping reads allocated by CSEM. Otherwise, please provide FASTQ format DNase-seq file.", sep = ""))
        file.bed <- dnaseFile
    }

  
    ## If it wasn't created before or given by the user, creates chrList from the bed file
    if(is.null(chrList)){
        system(paste("awk '{print $1}' ", file.bed, " | sort | uniq > ", outfileLoc, "/chrList_temp", sep=""))
        chrList <- as.character(read.table(paste(outfileLoc, "/chrList_temp", sep=""))$V1)
        system(paste("rm -rf ", outfileLoc, "/chrList_temp", sep=""))
    }
  
    # Filter out the alignments with posterior probability < AllocThres/1000
    system(paste("awk '$5 > ", AllocThres, " {print $0}' ", file.bed, ">", file.bed_AllocThres, sep=""))
  
    message( "Info: Reading the DNase-seq aligned read file and processing it into nucleotide-level files..." )
    .constructNeucleotideLevelCounts(outfileLoc, file.bed_AllocThres, fragL, outfileLoc, chrList, capping)
    countfile <- paste('_', file.bed_AllocThres, '_fragL', fragL, '_bin1.txt', sep='')
  
    message( "Info: Partitioning genome (DNase-seq) ..." )
    r <- .clusterPositions(countfile, chrList, outfileLoc, chrom.ref, outfile)
    dnaseKnots <- r[[1]]
    dnaseThres <- r[[2]]
    link <- paste(outfileLoc, '/', outfile, '_positions_cluster.txt', sep='')
    link_bychr <- vector('list', length(chrList))
    names(link_bychr) <- chrList
    for(i in chrList){
        link_bychr[[i]] <- paste(outfileLoc, '/', i, '_', outfile, '_positions_cluster.txt', sep='')
    }
    # Remove files if saveFiles=FALSE
    if(saveFiles=='FALSE' | saveFiles=='F'){
        system(paste('rm -rf ', file.sai, sep = ''))
        #system(paste("rm -rf ", file.bam, sep = ""))
        #system(paste("rm -rf ", file.sam, sep = ""))
        system(paste("rm -rf *sorted.bam*", sep = ""))
        #system(paste("rm -rf ", '*_', AllocThres, '.bed*', sep = ""))
        #system(paste("rm -rf *_output_Dnasequantile.txt*", sep = ""))
    }
  
    dnaseObject <- list(dnaseKnots = dnaseKnots, dnaseThres = dnaseThres, posLoc_bychr = link_bychr, chrList = chrList)
    return(dnaseObject)
}
