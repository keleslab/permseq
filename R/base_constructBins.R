# Read alignment files and construct bin-level files
.constructBins <- function( infile=NULL, fileFormat=NULL, outfileLoc="./", 
                           byChr=FALSE, useChrfile=FALSE, chrfile=NULL, excludeChr=NULL, 
                           fragLen=200, binSize=1, capping=0, perl = "perl" ){   
  # preprocessing perl script embedded in "mosaics/inst/Perl/"
  PET=FALSE
  if( PET == TRUE ) {
    script <- "process_readfiles_PET.pl"
    allFormat <- c( "eland_result", "sam" )
    allFormatName <- c( "Eland result", "SAM" )
  }else{
    script <- "process_readfiles_SET.pl"
    allFormat <- c( "eland_result", "eland_extended", "eland_export", 
                   "bowtie", "sam", "bed", "csem" )
    allFormatName <- c( "Eland result", "Eland extended", "Eland export", 
                       "Bowtie default", "SAM", "BED", "CSEM" )
  }
  #Check whether perl exists
  CMD <- paste( perl, "-v" )
  res <- system( CMD, intern = TRUE, ignore.stderr = TRUE )
  if( length(res) == 0 ){
    # cannot proceed if perl does not exist
    stop( "Perl is not found on your system! Either check $PATH if installed or please install Perl." )
  }else{
  # process read files into bin-level files if perl exists
  # check whether minimal options are missing      
    if( length(infile) != 1 || is.null(infile) )
      stop( "Please specify the name of the aligned read file!" )
          
    if( length(fileFormat) != 1 || is.null(fileFormat) )
      stop( "Please specify aligned read file format! Read '?constructBins' for supported file formats" )
    
      # check file format specification
    if( length(which(!is.na(match( fileFormat, allFormat )))) == 0 )
      stop( "Unsupported aligned read file format! Read '?constructBins' for supported file formats" )
      
      # if useChrfile is TRUE & excludeChr is NOT null, then ignore excludeChr  
    if( useChrfile & !is.null(excludeChr) ){
      message( "User set 'useChrfile' as TRUE and also provided 'excludeChr'." )
      message( "'excludeChr' argument will be ignored." )
      excludeChr <- NULL
    }
    # print out processing settings:
    # by default, set fragment length = 200, bin size = fragment length, capping = 0.        
    fileFormatName <- allFormatName[ match( fileFormat, allFormat ) ]
    cat( "------------------------------------------------------------\n" )
    cat( "Info: setting summary\n" )
    cat( "------------------------------------------------------------\n" )
    cat( "Name of aligned read file:", infile, "\n" )
    cat( "Aligned read file format:", fileFormatName, "\n" )
    cat( "Directory of processed bin-level files:", outfileLoc, "\n" )
    if ( useChrfile == TRUE ){ 
      cat( "Name of file for chromosome info: ", chrfile, "\n" )
    }
    if( !is.null(excludeChr) ){
      cat( "List of chromosomes to be excluded:", paste(excludeChr, collapse=", "), "\n" )
    }
    if( PET == FALSE ){
      cat( "Data type: Single-end tag (SET)\n" )
      cat( "Average fragment length:", fragLen, "\n" )
    }else{
      cat( "Data type: Paired-end tag (PET)\n" )
    }
    cat( "Bin size:", binSize, "\n" )
    if( capping > 0 ){
      cat( "Maximum number of reads allowed in each nucleotide:", capping, "\n" )
    }
    cat( "------------------------------------------------------------\n" )
    # get path to the perl code (unified script for all file formats)
    Fn.Path <- system.file( file.path("Perl",script), package="permseq")
    # process read file to bin-level files using perl codes
    message( "Info: reading the aligned read file and processing it into bin-level files..." )
    if( capping <= 0 ){
      capping <- 0
    }
    if( is.null(excludeChr) ){
      excludeChrVec <- ""
    }else{
      excludeChrVec <- paste( excludeChr, collapse=" " )
    }
    if( PET == TRUE ){
      CMD <- paste( perl, 
                   " ", Fn.Path,
                   " ", infile, 
                   " ", outfileLoc, 
                   " ", fileFormat,
                   " ", binSize, 
                   " ", capping,
                   " ", ifelse(byChr,"Y","N"),
                   " ", ifelse(useChrfile,"Y","N"),
                   " ", ifelse(!is.null(chrfile),chrfile,"-"),
                   " ", paste(excludeChrVec,collpase=" "), sep="" )
    }else{    
      CMD <- paste( perl, 
                   " ", Fn.Path,
                   " ", infile, 
                   " ", outfileLoc, 
                   " ", fileFormat, 
                   " ", fragLen, 
                   " ", binSize, 
                   " ", capping, 
                   " ", ifelse(byChr,"Y","N"),
                   " ", ifelse(useChrfile,"Y","N"),
                   " ", ifelse(!is.null(chrfile),chrfile,"-"),
                   " ", paste(excludeChrVec,collpase=" "), sep="" )
    }
    res <- system( CMD, intern = TRUE )
    message( "Info: done!" )   
    # print out processing results
    infilename <- basename( infile )
    # extract only filename from infile
    byChr = FALSE
    if( PET == TRUE ){    
      if( byChr ){
        outfileName <- list.files( path = outfileLoc, pattern = paste(infilename, "_bin", binSize, "_.*.txt", sep="") )
      }else{
        outfileName <- paste(infilename, "_bin", binSize, ".txt", sep="")
      }
    }else{
      if( byChr ){
        outfileName <- list.files( path = outfileLoc, pattern = paste(infilename, "_fragL", fragLen, "_bin", binSize, "_.*.txt", sep="") )
      }else{
        outfileName <- paste(infilename, "_fragL", fragLen, "_bin", binSize, ".txt", sep = "")
      }
    }
    
    cat( "------------------------------------------------------------\n" )
    cat( "Info: processing summary\n" )
    cat( "------------------------------------------------------------\n" )    
    cat( "Directory of processed bin-level files:", outfileLoc, "\n" )
    if( byChr ){
      cat( "List of processed bin-level files:\n" )
      for( i in 1:length(outfileName) ){                
        cat( paste(i, ". ", outfileName[i], "\n", sep="") )
      }
    }else{
      cat( "Processed bin-level file: ", outfileName, "\n", sep="" )   
    }
    cat( "------------------------------------------------------------\n" )
  }
}
