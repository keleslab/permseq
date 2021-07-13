
setMethod(f = "[", signature = "Prior", definition=function(x, i, drop="missing"){
            if(i == "dnaseName"){return(x@dnaseName)}else{}
      
            if(i == "dnaseAlign"){return(x@dnaseAlign)}else{}
            if(i == "dnaseKnots"){return(x@dnaseKnots)}else{}
            if(i == "dnaseThres"){return(x@dnaseThres)}else{}
            if(i == "posLoc_bychr"){return(x@posLoc_bychr)}else{}
            if(i == "dnaseHistone"){return(x@dnaseHistone)}else{}
            
            if(i == "histoneName"){return(x@histoneName)}else{}
            if(i == "histoneNum"){return(x@histoneNum)}else{}
            if(i == "histoneAlign"){return(x@histoneAlign)}else{}
            if(i == "histoneGrpL"){return(x@histoneGrpL)}else{}

            if(i == "chipName"){return(x@chipName)}else{}
            if(i == "chipNum"){return(x@chipNum)}else{}
            if(i == "chipAlign"){return(x@chipAlign)}
            if(i == "chipSAM"){return(x@chipSAM)}else{}
            if(i == "chipAllocate"){return(x@chipAllocate)}else{}
            if(i == "chipUni"){return(x@chipUni)}else{}
            if(i == "chipAlignFormat"){return(x@chipAlignFormat)}else{}


            if(i == "dataNum"){return(x@dataNum)}else{}
            if(i == "fragL"){return(x@fragL)}else{}
            if(i == "chrList"){return(x@chrList)}else{}
            if(i == "bowtieInfo"){return(x@bowtieInfo)}else{}
            if(i == "bwaInfo"){return(x@bwaInfo)}else{}
            if(i == "csemDir"){return(x@csemDir)}else{}
            if(i == "picardDir"){return(x@picardDir)}else{}
            if(i == "outfileLoc"){return(x@outfileLoc)}else{}
            if(i == "prior"){return(x@prior)}else{}
            if(i == "chrom.ref"){return(x@chrom.ref)}else{}
            


          }
          )



##' Extract part of the Prior object("[[")
##' 
##' @param x Prior object
##' @param i Name for the extraction part
##' @param drop Missing if i is absent
##' @rdname subsubpriorMethods
setMethod(f = "[[", signature = "Prior", definition=function(x, i, drop="missing"){

  if(i == "dnaseName"){return(x@dnaseName)}else{}

  if(i == "dnaseAlign"){return(x@dnaseAlign)}else{}
  if(i == "dnaseKnots"){return(x@dnaseKnots)}else{}
  if(i == "dnaseThres"){return(x@dnaseThres)}else{}
  if(i == "posLoc_bychr"){return(x@posLoc_bychr)}else{}
  if(i == "dnaseHistone"){return(x@dnaseHistone)}else{}
  
  if(i == "histoneName"){return(x@histoneName)}else{}
  if(i == "histoneNum"){return(x@histoneNum)}else{}
  if(i == "histoneAlign"){return(x@histoneAlign)}else{}
  if(i == "histoneGrpL"){return(x@histoneGrpL)}else{}

  if(i == "chipName"){return(x@chipName)}else{}
  if(i == "chipNum"){return(x@chipNum)}else{}
  if(i == "chipAlign"){return(x@chipAlign)}
  if(i == "chipSAM"){return(x@chipSAM)}else{}
  if(i == "chipAllocate"){return(x@chipAllocate)}else{}
  if(i == "chipUni"){return(x@chipUni)}else{}
  if(i == "chipAlignFormat"){return(x@chipAlignFormat)}else{}

  
  if(i == "dataNum"){return(x@dataNum)}else{}
  if(i == "fragL"){return(x@fragL)}else{}
  if(i == "chrList"){return(x@chrList)}else{}
  if(i == "bwaInfo"){return(x@bwaInfo)}else{}
  if(i == "bowtieInfo"){return(x@bowtieInfo)}else{}
  if(i == "csemDir"){return(x@csemDir)}else{}
  if(i == "picardDir"){return(x@picardDir)}else{}
  if(i == "outfileLoc"){return(x@outfileLoc)}else{}
  if(i == "prior"){return(x@prior)}else{}
  if(i == "chrom.ref"){return(x@chrom.ref)}else{}
          
            
}
          )

##' Print function for Prior class
##'
##' print all the information within Prior class
##' @title Print
##' @param x Prior object
##' @return Screen print out
##' @rdname priorMethods
setMethod("print", "Prior",
          function(x){
            cat("*******************************************************\n")
            cat("********Print all the information in 'Prior' Object****\n")
            cat("*******************************************************\n")
            cat("** DNase-seq Related Information:\n")
      
            cat("Dataset used as DNase-seq:\n")
            print(x@dnaseName)
            cat("\n")
            cat("DNase-seq alignment information:\n")
            print(x@dnaseAlign)
            cat("\n")
            cat("B-spline function knots:\n")
            print(x@dnaseKnots)
            cat("\n")
            cat("DNase-seq group (showing the first 30):\n")
            print(x@dnaseThres[1:min(30, length(x@dnaseThres))])
            cat("\n")
            cat("DNase-seq candidates when on DNase-seq dataset available:\n")
            print(x@dnaseHistone)
            cat("\n")
            cat("Position of files containing which segments of genome are in which group:\n")
            print(x@posLoc_bychr)
            cat("\n")
            
            cat("*****************************************************\n")
            cat("** Histone Related Information:\n")
            cat("Histone dataset(s) name:\n")
            print(x@histoneName)
            cat("\n")
            cat("Number of histone dataset(s):\n")
            print(x@histoneNum)
            cat("\n")
            cat("Histone dataset(s) alignment information:\n")
            print(x@histoneAlign)
            cat("\n")
            cat("Selected histone dataset(s) after group lasso variable selection:\n")
            print(x@histoneGrpL)
            cat("\n")
            
            cat("*****************************************************\n")
            cat("** ChIP-seq Related Information:\n")
            cat("ChIP-seq dataset(s) name:\n")
            print(x@chipName)
            cat("\n")
            cat("Number of ChIP-seq dataset(s):\n")
            print(x@chipNum)
            cat("\n")
            cat("ChIP-seq alignment information:\n")
            print(x@chipAlign)
            cat("\n")
            cat("ChIP-seq aligned by bowtie in SAM format:\n")
            print(x@chipSAM)
            cat("\n")
            cat("ChIP-seq uni-reads alignment in BED format:\n")
            print(x@chipUni)
            cat("\n")
            cat("ChIP-seq aligned by Permseq in BAM format:\n")
            print(x@chipAllocate)
            cat("\n")
            cat("ChIP-seq aligned by Permseq in other format:\n")
            print(x@chipAlignFormat)
            cat("\n")
            
            cat("*****************************************************\n")
            cat("** Other Parameter and Directory Information:\n")
            cat("Total of dataset(s):\n")
            print(x@dataNum)
            cat("\n")
            cat("Fragment length:\n")
            print(x@fragL)
            cat("\n")
            cat("Bowtie parameter and path related information:\n")
            print(x@bowtieInfo)
            cat("\n")
            cat("CSEM directory:\n")
            print(x@csemDir)
            cat("\n")
            cat("Picard directory:\n")
            print(x@picardDir)
            cat("\n")
            cat("Directory to store the output files:\n")
            print(x@outfileLoc)
            cat("\n")
            cat("Location where the prior is saved:\n")
            print(x@prior)
            cat("\n")
            cat("Chromesome summary file:\n")
            print(x@chrom.ref)
            cat("\n")
              

            cat("********************************************************\n")
            cat("****End Print of all the information in 'Prior' object**\n")
 
            cat("********************************************************\n")           
          })


##' Show function for Prior class
##'
##' show the selected important information in Prior object
##' @title Show
##' @param object Prior object 
##' @return Screen print out
##' @rdname priorMethods
setMethod("show", signature(object = "Prior"),
          definition = function(object){
            

            cat("*****************************************************\n")
            cat("**********Prior Object Minimal Information***********\n")
            cat("*****************************************************\n")
            cat("\n")
            cat("****DNase-seq Minimal Information:\n")
            cat("\n")

            cat("** Name of dataset used as DNase-seq data:\n")
            if(length(object@dnaseName) == 0){
              print("NULL")
            }else{
              print(object@dnaseName)
            }
            cat("\n")

            cat("** 90, 99, 99.9th percentile of DNase-seq read count:\n")
            if(length(object@dnaseKnots) == 0){
              print("NULL")
            }else{
            
              print(object@dnaseKnots)
            }
            cat("\n")


            cat("** DNase threshold matrix:\n")
            if(length(object@dnaseThres) == 0){

              print("NULL")
            }else{
              cat("Head 10 of dnaseThres:\n")
              print(object@dnaseThres[1:10])
              cat('......\n')
              cat("Tail 10 of dnaseThres:\n ")
              print(object@dnaseThres[(length(object@dnaseThres) - 10) : length(object@dnaseThres)])
            }
            cat("\n")
            

            cat("** In OnlyHistone, the Potential DNase aligned by Histone data:\n")
            if(length(object@dnaseHistone) == 0){
              print("NULL")
            }else{
              
              potentialDnase <- object@dnaseHistone
              potentialNames <- names(object@dnaseHistone)
              for(i in 1:object@dataNum){
                potentialDnase[[potentialNames[i]]]$dnaseThres <- potentialDnase[[potentialNames[i]]]$dnaseThres[1:min(100, length(potentialDnase[[potentialNames[i]]]$dnaseThres))]
              }
              print(potentialDnase)
              
            }
            cat("\n")

            cat("*****************************************************\n")
            cat("****Histone Minimal Information:\n")
            cat("\n")

            cat("** Histone Data Name:\n")
            if(length(object@histoneName) == 0){
              print("NULL")
            }else{
              print(object@histoneName)
            }
            cat("\n")
           
            cat("** Selected Histone Dataset(s):\n")
            if(length(object@histoneGrpL) == 0){
              print("NULL")
            }else{
            
              print(object@histoneGrpL)
            }
            cat("\n")


            cat("*****************************************************\n")
            cat("****ChIP-seq Minimal Information:\n")
            cat("\n")
            cat("** ChIP-seq Name:\n")
            if(length(object@chipName) == 0){
              print("NULL")
            }else{
              print(object@chipName)
            }
            cat("\n")
           
            
            cat("** ChIP-seq Output Directory and Format:\n")
            if(length(object@chipAlignFormat) == 0){
              print("NULL")
            }else{
            
              print(object@chipAlignFormat)
            }
            cat("\n")

            cat("*****************************************************\n")
            cat("****Other Information:\n")

            cat("\n")

            cat("** Total number of DNase-seq and Histone datasets:\n")
            if(length(object@dataNum) == 0){
              print("NULL")
            }else{
            
              print(object@dataNum)
            }

            cat("\n")

            cat("** Directory where prior is saved:\n")
            if(length(object@prior) == 0){
              print("NULL")
            }else{
            
              print(object@prior)
            }
            cat("\n")

            cat("** Directory where the chromosome reference information is saved\n")
            if(length(object@chrom.ref) == 0){
              print("NULL")
            }else{
            
              print(object@chrom.ref)
            }
            cat("\n")

            
            cat("*****************************************************\n")
            cat("*************End of the Core Information*************\n")
            cat("*****************************************************\n")
            
            
          })



##' Plot the histone versus ChIP read counts 
##'
##' plot the histone versus ChIP read counts
##' @title Plot
##' @param x Prior object
##' @param y Other parameter
##' @param ... other parameter
##' @return Histone versus ChIP read counts 
##' @rdname priorMethods
setMethod(f="plot", signature = "Prior",  definition = function(x, y, ...){
  
  if(length(x@chipSAM) == 0){
    stop("ChIP-seq alignment file(s) are not available now. Please plot after priorGenerate process!")
    
  }else{
    if(length(x@posLoc_bychr) == 0){
      stop("DNase-seq read count is not available now. Please provide the DNase-seq data or wait until histone dataset is selected to work as DNase-seq!")
    }else{
      
      
      chipSAM <- x@chipSAM
      outfileLoc <- x@outfileLoc
      chipFileName <- gsub(".*/(.*).sam", "\\1", chipSAM)
      outfile_chip <- paste(chipFileName, ".sam", sep = "")
      dnaseThres <- x@dnaseThres
      dnaseKnots <- x@dnaseKnots
      
      outfile_chipmean <- paste("chip", 1:length(outfile_chip), "_chipmean", sep="")

      posLoc_bychr <- vector('list', length(x@dnaseName))
      names(posLoc_bychr) <- x@dnaseName
      posLoc_bychr[[length(x@dnaseName)]] <- x@posLoc_bychr
      
      #calculate averaged chip read counts according to object (dnase or histone information)
      if(!file.exists((paste(outfileLoc, '/', names(posLoc_bychr)[1], '_', outfile_chipmean[1], sep='')))){
        .chipMeanCounts(x, posLoc_bychr, chipSAM, outfile_chip, outfileLoc, outfile_chipmean)
      }
     #results are reported in outfile_chipmean
      
####################figure generate######################
      
      .fitPlot2 = function(reps, dnaseThres, dnaseKnots, name, ylim, ...)
        {
          par(mar = c(5, 5, 4, 2) + 0.1)
          
          for (i in 1:length(reps)) {
          
            a <- reps[[i]][2, ]/reps[[i]][1, ]
            matplot(dnaseThres, t(a), pch = 20, cex = 0.5, ylab = "Average ChIP read count", xlab = "DNase read count", main=name, ylim=ylim, ...)
            
            for (j in 1:length(dnaseKnots))
              abline(v = dnaseKnots[j], col = "blue", lwd = 2)
              #legend(x = "top", legend = c("Observed data", "Fitted value"),
              #   col = c("black", "red"), lwd = c(1, 3), lty = c(0,
              #      2), pch = c(1, NA), bty = "n")
          }
        }
      
      #read in mean chip counts accoding to each DNase or histone read counts
      ylim <- vector('list', length(outfile_chipmean))
      for(i in 1:length(outfile_chipmean))
        ylim[[i]] <- c(0,0)
      reps <- vector('list', length(posLoc_bychr))
      for(i in 1:length(posLoc_bychr)){
        reps[[i]] <- vector('list', length(outfile_chipmean))
        for(j in 1:length(outfile_chipmean)){
          reps[[i]][[j]] <- read.table(paste(outfileLoc, '/', names(posLoc_bychr)[i], '_', outfile_chipmean[j], sep=''))  #read in averaged Chip read counts
          ylim[[j]][1] <- min(c(ylim[[j]][1], min(unlist(reps[[i]][[j]][2,]/reps[[i]][[j]][1,]))))
          ylim[[j]][2] <- max(c(ylim[[j]][2], max(unlist(reps[[i]][[j]][2,]/reps[[i]][[j]][1,]))))
          
                                        #find the min and max average chip counts
          
        }
      }
####################plot ##################################
                                        #pdf(paste(outfileLoc, 'marginal_plot_test.pdf', sep = ""))
      par(mfrow=c(length(outfile_chipmean), length(posLoc_bychr)))
      for(j in 1:length(outfile_chipmean)){
        for(i in 1:length(posLoc_bychr)){
          .fitPlot2(list(reps[[i]][[j]]), dnaseThres, dnaseKnots,name= paste(names(posLoc_bychr)[i], "_", chipFileName[j], sep = ""), ylim=ylim[[j]], ...)
          

        }
        
        
      }
     
    }
  }
}
          )
          
          


##'
##' Set value to the extract part of the Prior object
##' 
##' @param x Prior object
##' @param i Name for the extraction part
##' @param value The updated value
##' @rdname subpriorMethods
setReplaceMethod(f = "[", signature = "Prior", definition=function(x,i,value){
  


            if(i == "dnaseName"){x@dnaseName <- value}else{}
          
            if(i == "dnaseAlign"){x@dnaseAlign <- value}else{}
            if(i == "dnaseKnots"){x@dnaseKnots <- value}else{}
            if(i == "dnaseThres"){x@dnaseThres <- value}else{}
            if(i == "dnaseHistone"){x@dnaseHistone <- value}else{}

            if(i == "histoneName"){x@histoneName <- value}else{}
            if(i == "histoneNum"){x@histoneNum <- value}else{}
            if(i == "histoneAlign"){x@histoneAlign <- value}else{}
            if(i == "histoneGrpL"){x@histoneGrpL <- value}else{}

            if(i == "chipName"){x@chipName <- value}else{}
            if(i == "chipNum"){x@chipNum <- value}else{}
            if(i == "chipAlign"){x@chipAlign <- value}else{}
            if(i == "chipSAM"){x@chipSAM <- value}else{}
            if(i == "chipAllocate"){x@chipAllocate <- value}else{}
            if(i == "chipUni"){x@chipUni <- value}else{}
            if(i == "chipAlignFormat"){x@chipAlignFormat <- value}else{}


                 
            if(i == "dataNum"){x@dataNum <- value}else{}
            if(i == "fragL"){x@fragL <- value}else{}
            if(i == "chrList"){x@chrList <- value}else{}
            if(i == "bowtieInfo"){x@bowtieInfo <- value}else{}
            if(i == "bwaInfo"){x@bwaInfo <- value}else{}
            if(i == "csemDir"){x@csemDir <- value}else{}
            if(i == "picardDir"){x@picardDir <- value}else{}
            if(i == "outfileLoc"){x@outfileLoc <- value}else{}
            if(i == "prior"){x@prior <- value}else{}
            if(i == "posLoc_bychr"){x@posLoc_bychr <- value}else{}
            if(i == "chrom.ref"){x@chrom.ref <- value}else{}
            return(x)

})

##'
##' Set value to the extract part of the Prior object("[[")
##' 
##' @param x Prior object
##' @param i Name for the extraction part
##' @param value The updated value
##' @rdname subsubpriorMethods
setReplaceMethod(f = "[[", signature = "Prior", definition=function(x,i,value){
  

            if(i == "dnaseName"){x@dnaseName <- value}else{}
  
            if(i == "dnaseAlign"){x@dnaseAlign <- value}else{}
            if(i == "dnaseKnots"){x@dnaseKnots <- value}else{}
            if(i == "dnaseThres"){x@dnaseThres <- value}else{}
            if(i == "dnaseHistone"){x@dnaseHistone <- value}else{}
            
            if(i == "histoneName"){x@histoneName <- value}else{}
            if(i == "histoneNum"){x@histoneNum <- value}else{}
            if(i == "histoneAlign"){x@histoneAlign <- value}else{}
            if(i == "histoneGrpL"){x@histoneGrpL <- value}else{}

            if(i == "chipName"){x@chipName <- value}else{}
            if(i == "chipNum"){x@chipNum <- value}else{}
            if(i == "chipAlign"){x@chipAlign <- value}else{}
            if(i == "chipSAM"){x@chipSAM <- value}else{}
            if(i == "chipAllocate"){x@chipAllocate <- value}else{}
            if(i == "chipUni"){x@chipUni <- value}else{}
            if(i == "chipAlignFormat"){x@chipAlignFormat <- value}else{}


                 
            if(i == "dataNum"){x@dataNum <- value}else{}
            if(i == "fragL"){x@fragL <- value}else{}
            if(i == "chrList"){x@chrList <- value}else{}
            if(i == "bowtieInfo"){x@bowtieInfo <- value}else{}
            if(i == "bwaInfo"){x@bwaInfo <- value}else{}
            if(i == "csemDir"){x@csemDir <- value}else{}
            if(i == "picardDir"){x@picardDir <- value}else{}
            if(i == "outfileLoc"){x@outfileLoc <- value}else{}
            if(i == "prior"){x@prior <- value}else{}
            if(i == "posLoc_bychr"){x@posLoc_bychr <- value}else{}
            if(i == "chrom.ref"){x@chrom.ref <- value}else{}
            return(x)



})

##' return information about the individual slots in "Prior" object.
##'
##' names()=slotNames()
##' @title names
##' @param x Prior object
##' @return Names of the individual slots in "Prior" object
##' @rdname priorMethods
setMethod("names", "Prior",
	function(x){
	 return(slotNames(x))	
	}
          )


##'
##' return Bowtie Alignment Information
##' @title summary
##' @param object Prior object
##' @return Alignment details of DNase-seq, Histone, ChIP-seq
##' @rdname priorMethods
setMethod(f="summary", signature = "Prior",  definition = function(object){

  x = object
  if(length(x@dnaseAlign)!=0){
    AlignInfo <- x@dnaseAlign
    cat("***********************************************************************************\n")
    cat(paste("Alignment Summary for DNase-seq Datasets ", x@dnaseName, ":\n", sep = ""))
    cat("***********************************************************************************\n")
    table <- matrix(c(AlignInfo$readNum, AlignInfo$readAlignedNum, AlignInfo$readAlignedRate, AlignInfo$readFailedNum, AlignInfo$readFailedRate, AlignInfo$AlignNum), ncol = 1)
    colnames(table) <- c("Alignment Information")
    rownames(table) <- c("Number of reads processed:", "Number of reads with at least one reported alignment:", "Percentage of reads with at least one reported alignment(%):", "Number of reads failed to align:", "Percentage of reads failed to align(%):", "Total number of alignment reported:")
    
    print(table) 
     

  }
  
if(length(x@histoneAlign) != 0){
  for(i in 1:length(x@histoneName)){
    if(length(x@histoneAlign[[i]])!=0){
                
      AlignInfo <- x@histoneAlign[[x@histoneName[i]]]
      cat("***********************************************************************************\n")
      cat(paste("Alignment Summary for Histone Dataset ", x@histoneName[i], " :\n", sep = ""))
      
      cat("***********************************************************************************\n")
      table <- matrix(c(AlignInfo$readNum, AlignInfo$readAlignedNum, AlignInfo$readAlignedRate, AlignInfo$readFailedNum, AlignInfo$readFailedRate, AlignInfo$AlignNum), ncol = 1)
      colnames(table) <- c("Alignment Information")
      rownames(table) <- c("Number of reads processed:", "Number of reads with at least one reported alignment:", "Percentage of reads with at least one reported alignment(%):", "Number of reads failed to align:", "Percentage of reads failed to align(%):", "Total number of alignment reported:")
    
      print(table)
    
                
    }
    
    
  }
}
  
 if(length(x@chipAlign) != 0){             
  for(i in 1:length(x@chipName)){
    if(length(x@chipAlign[[i]])!=0){
       
      AlignInfo <- x@chipAlign[[x@chipName[i]]]
      cat("***********************************************************************************\n")
      cat(paste("Alignment Summary for ChIP-seq Dataset ", x@chipName[i],":", sep = ""))
      cat("\n")      
      cat("***********************************************************************************\n")
      table <- matrix(c(AlignInfo$readNum, AlignInfo$readAlignedNum, AlignInfo$readAlignedRate, AlignInfo$readFailedNum, AlignInfo$readFailedRate, AlignInfo$AlignNum), ncol = 1)
      colnames(table) <- c("Alignment Information")
      rownames(table) <- c("Number of reads processed:", "Number of reads with at least one reported alignment:", "Percentage of reads with at least one reported alignment(%):", "Number of reads failed to align:", "Percentage of reads failed to align(%):", "Total number of alignment reported:")
    
      print(table)

      
      
 

                
    }
    
    
  }
}
  
})
