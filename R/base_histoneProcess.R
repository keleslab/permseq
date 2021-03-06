.histoneProcess = function(histoneFile, fragL, AllocThres, chrList, chrom.ref, capping, outfileLoc="./", outfile, bowtieIndex, csemDir, bowtieDir, vBowtie, mBowtie, pBowtie, bwaDir, bwaIndex, nBWA, oBWA, tBWA, mBWA, picardDir, saveFiles){

    if(!file.exists(outfileLoc)){
        system(paste('mkdir ',outfileLoc,sep=''))
    }
    
    setwd(outfileLoc)
    
    histoneFormat <- tolower(sub(".*\\.", "", histoneFile))
                                        # Define filenames and extensions 
    filename <- strsplit(histoneFile, "/")[[1]]
    filename <- filename[length(filename)]
    filename <- strsplit(filename, paste('.', sub(".*\\.", "", filename), sep=''))[[1]][1]
    file.sai <- paste(filename, '.sai', sep='')
    file.sam <- paste(filename, '.sam', sep='')
    file.bam <- paste(filename, '.bam', sep='')
    file.bed <- paste(filename, '.bed', sep='')
    file.bed_AllocThres <- paste(filename, '_', AllocThres, '.bed', sep='')
  
    ## Deal with different Histone file formats (fastq, sam, bam or bed)
    if(histoneFormat != "bed"){
        if(histoneFormat != "bam"){
            if(histoneFormat != "sam"){#Histone file is fastq format - will go through: alignment - CSEM multi-reads allocation steps
                if(!is.null(bowtieIndex)){
                    message( "Info: Aligning Histone ChIP-seq reads using Bowtie." )
                    ## Generates SAM file using Bowtie.
                    system(paste(bowtieDir,'/bowtie -q -v ',vBowtie,' -a -m ',mBowtie,' -p ',pBowtie,' --sam ', bowtieIndex,' ', histoneFile,' ',file.sam, " 2>&1 | tee ", outfileLoc, "/priorProcessHistoneBowtie_temp.txt", sep=''))  
                    
                }else{
                    message( "Info: Aligning Histone ChIP-seq reads using BWA." )
                    ## Generates the SAM file using BWA.
                    system(paste(bwaDir, '/bwa aln -n ', nBWA, ' -o ', oBWA, ' -t ', tBWA, ' ', bwaIndex, ' ', histoneFile, ' >', file.sai, sep=''))
                    system(paste(bwaDir, '/bwa samse -n ', mBWA, ' ', bwaIndex, ' ', file.sai, ' ', histoneFile, ' >', file.sam, '.multiOneLine', sep=''))

                    system(cat("awk '{split($0, tag,\"XA\");split($1, head, \"\"); if (head[1] == \"@\") print $0; else if ($5 >24) print $0; else if (tag[2] != \"\") print $0;}' ", file.sam, ".multiOneLine | ", bwaDir, "/xa2multi.pl >", file.sam, ".tmp\n", sep = ""))
                     
                    system(cat("cat ", file.sam, ".tmp | awk 'BEGIN {FS=\"\\t\" ; OFS=\"\\t\"} ! /^@/ && $6!=\"*\" { cigar=$6; gsub(\"[0-9]+D\",\"\",cigar); n = split(cigar,vals,\"[A-Z]\"); s = 0; for (i=1;i<=n;i++) s=s+vals[i]; seqlen=length($10) ; if (s!=seqlen) {print $0} }' >", file.sam, ".tmp.badCIGAR\n", sep = ""))

                    system(cat("if [ $(cat ", file.sam, ".tmp.badCIGAR | wc -l) -gt 0 ]; then\ncat ", file.sam, ".tmp  | grep -v -F -f ", file.sam, ".tmp.badCIGAR >", file.sam, "\nelse\nmv ", file.sam, ".tmp ", file.sam, "\nfi\n", sep = "")) 


                    ## system(paste(bwaDir, '/bwa samse -n ', mBWA, ' ', bwaIndex, ' ', file.sai, ' ', histoneFile, ' | ', bwaDir, '/xa2multi.pl >', file.sam, sep=''))
                }
            }else{#Histone file is in SAM format -  can skip alignment step and directly allocate multi-reads using CSEM
        
                print(paste("Histone file ", filename, " is in SAM format. The preprocessed aligned SAM file should contain multi-mapping reads. Otherwise, please provide FASTQ format Histone file.", sep = ""))
                file.sam <- histoneFile
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
            

            message( "Info: Allocating Histone ChIP-seq multi-reads sam file using CSEM." )
            if(!is.null(bowtieIndex)){

                system(paste(csemDir, "/run-csem --sam -p ", pBowtie, ' ', file.sam, ' ', fragL, ' ', filename, sep=''))
                
            }else{
                
                system(paste(csemDir, "/run-csem --sam -p ", tBWA, ' ', file.sam, ' ', fragL, ' ', filename, sep=''))

            }
      
        }else{#Histone file is in BAM format - It is better to be the CSEM alignemnt results

            print(paste("Histone file ", filename, " is in SAM format. The preprocessed aligned BAM file should contain multi-mapping reads allocated by CSEM. Otherwise, please provide FASTQ format Histone file.", sep = ""))
            file.bam <- histoneFile
        }
    
        message( "Info: Transfering Histone ChIP-seq multi-reads alignment bam file into bed file using CSEM" )
        system(paste(csemDir,"/csem-generate-input --bed ",file.bam,' ',filename,sep='')) 
        
    }else{#Histone ChIP-seq file is in BED format - It is better to be the CSEM alignment results and will be used directly for prior construction

        print(paste("Histone file ", filename, " is in SAM format. The preprocessed aligned BED file should contain multi-mapping reads allocated by CSEM. Otherwise, please provide FASTQ format Histone file.", sep = ""))
        file.bed <- histoneFile

    }
  
  # Filter out the alignments with posterior probability < AllocThres/1000
    system(paste("awk '$5 > ", AllocThres, " {print $0}' ", file.bed, ">", file.bed_AllocThres, sep=""))
  
    message( "Info: Reading the Histone ChIP-seq aligned read file and processing it into nucleotide-level files..." )
    .constructNeucleotideLevelCounts(outfileLoc, file.bed_AllocThres, fragL, outfileLoc, chrList, capping)
    countfile <- paste('_', file.bed_AllocThres, '_fragL', fragL, '_bin1.txt', sep='')
    
    message( "Info: Partitioning genome (Histone) ..." )
    r <- .clusterPositions_3cluster(countfile, chrList, outfileLoc, chrom.ref, outfile)
    link <- vector('list', length(chrList))
    names(link) <- chrList
    for(i in chrList){
        link[[i]] <- paste(outfileLoc, '/', i, '_', outfile, '_positions_3cluster.txt', sep='')
    }
    ## Remove files if saveFiles=FALSE
    if(saveFiles == 'FALSE' | saveFiles == 'F'){
        system(paste('rm -rf ', file.sai, sep = ''))
        ##system(paste("rm -rf ", file.bam, sep = ""))
        ##system(paste("rm -rf ", file.sam, sep = ""))
        system(paste("rm -rf *sorted.bam*", sep = ""))
        ##system(paste("rm -rf ", '*_', AllocThres, '.bed*', sep = ""))
        ##system(paste("rm -rf *_output_Dnasequantile.txt*", sep = ""))
    }
    return(link)
}
