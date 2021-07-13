\name{permseq-package}
\alias{permseq}
\docType{package}
\title{Mapping protein-DNA interactions in highly repetitive regions of genomes with prior-enhanced read mapping.}
\description{
  Allocating multi-reads using the Multinomial-Dirichlet model.

  This package contains five main functions,\cr
  \code{\link{priorProcess}} \cr
  \code{\link{priorHistone_init}}\cr
  \code{\link{priorHistone_multi}}\cr
  \code{\link{priorGenerate}}\cr
  \code{\link{readAllocate}}\cr

  
  \code{\link{priorProcess}} constructs DNase and/or Histone cell-line specific module, it first constructs  nucleotide-level count files (DNase) and trinary indicators (Histone) from the aligned read file and constructs an object to be used in  \code{priorGenerate}.

  
  \code{\link{priorHistone_init}} processes Histone data  and select one Histone to be used as DNase data, when DNase dataset is not avaiable.
  
  \code{\link{priorHistone_multi}} processes Histone  data and generate module for further prior construction. 


  \code{\link{priorGenerate}} calculates averaged ChIP-seq read counts at different DNase and/or Histone counts group  and generates priors for ChIP data, it imports the object generated from \code{priorProcess} or  \code{priorHistone_multi}.

  \code{\link{readAllocate}} will allocate ChIP-seq multi-reads using prior information imported from the  object generated  by \code{priorGenerate}.

}
\author{
Xin Zeng, Sunduz Keles, Maria Constanza Rojo, Ye Zheng.\cr


Maintainer: Ye Zheng <yezheng@stat.wisc.edu>, M. Constanza Rojo-Alfaro <rojo@stat.wisc.edu>

}
\references{
Mapping protein-DNA interactions in segmental duplication and highly
repetitive regions of the genomes with prior-enhanced read mapping.
}
\examples{

\dontrun{
#When DNase-seq is available:
priorProcess_DNase <- priorProcess(dnaseFile="DNase_partial.fastq",
histoneFile=NULL, dnaseName="dnase", histoneName=NULL, fragL=200,
AllocThres=900, chrList=c("chr22"), capping=0, outfileLoc="./",
outfile="only_dnase", bowtieDir="./bowtie-1.1.1", bowtieIndex="./hg19",
vBowtie=2, mBowtie=99, pBowtie=8, csemDir="./csem-2.3",
picardDir=./picard.jar, saveFiles=TRUE)

#When DNase-seq and Histone are both available:
priorProcess_DNaseHistone <- priorProcess(dnaseFile="DNase_partial.fastq",
histoneFile=c("h2az_partial.fastq", "h3k4me1_partial.fastq"),
dnaseName="dnase", histoneName=c("h2az", "h3k4me1"), fragL=200,
AllocThres=900, chrList=c("chr22_partial"), capping=0,
outfileLoc="./", outfile="dnase_hiistone", bowtieIndex="hg19",
csemDir="./csem-2.3", picardDir=./picard.jar, bowtieDir="./bowtie-1.1.1", vBowtie=2,
mBowtie=99, pBowtie=8, saveFiles=TRUE)

#When DNase data not available:
priorHistone_init <- priorHistone_initial(histoneFile=c("h2az.fastq",
"h3k4me1.fastq"), histoneName=c("h2az", "h3k4me1"), fragL=200,
binSize=200, AllocThres=900, chrList=c("chr19", "chr22"), capping=0, o
utfileLoc="./", chipFile=c("chip_rep1.fastq", "chip_rep2.fastq"),
bowtieIndex="./hg19", csemDir="./csem-2.3", picardDir=./picard.jar, bowtieDir="./bowtie-1.1.1",
vBowtie=2, mBowtie=99, pBowtie=8, saveFiles=TRUE)

#select one data from the marginal plot
priorHistone_multi <- priorHistone_multi(priorHistone, dnaseIndex=1,
outfileLoc="./", outfile="twodata_histone")

#Generate Prior
priorGenerate <- priorGenerate(priorProcess_object,
chipFile="chip_partial.fastq", maxNumHistone=6, outfileLoc="./")

#Allocate the read
result <- readAllocate(priorGenerate, outfileLoc="./",
outputFormat="tagAlign", chipThres=500, chipFile=NULL, bowtieIndex=NULL,
csemDir=NULL, bowtieDir=NULL, vBowtie=NULL, mBowtie=NULL, pBowtie=NULL,
fragL=NULL)
}
}
