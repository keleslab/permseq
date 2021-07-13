\docType{class}
\name{priorClass}
\alias{priorClass}
\title{An S4 "Prior" class containing prior information for the ChIP-seq alignment.}
\arguments{

\item{dnaseName}{Names of datasets that are used as DNase-seq,
  especially in the "only Histone" situation, dnaseName is the selected
  histone ChIP-seq dataset.}
\item{dnaseAlign}{DNase-seq alignment summary information from bowtie.}
\item{dnaseKnots}{A vector of knots points for the B-spline
  functions. They are the 90, 99 and 99.9th percentiles of read counts.}
\item{dnaseThres}{A vector of DNase-seq group created to generate aggregated ChIP
data. After alignment, positions which have the same DNase-seq read count are clustered into
  one group. \code{dnaseThres} is the corresponding read count number
  in each group.}
\item{posLoc_bychr}{Location of the files containing the group index of
  each segment of the genome.}
\item{dnaseHistone}{The input histone ChIP-seq  files which are the candidates for the DNase-seq while running priorProcess procedure.}
\item{histoneName}{Name of histone ChIP-seq dataset(s). If no giving values, histoneName is set as a vector of index number 
\code{1:length(histoneFile)}.}
\item{histoneNum}{Number of histone ChIP-seq dataset(s).}
\item{histoneAlign}{Histone alignment summary information from bowtie.}
\item{histoneGrpL}{Histone dataset(s) selected to be included in the model after variable selection.}
\item{chipName}{Name of ChIP-seq dataset(s).}
\item{chipNum}{Number of ChIP-seq dataset(s).}
\item{chipAlign}{ChIP-seq alignment summary information from bowtie.}
\item{chipSAM}{Location of aligned ChIP-seq in SAM format.}
\item{chipAllocate}{Location of the final allocated ChIP-seq output
  by Permseq in BAM format.}
\item{chipUni}{Location of the aligned ChIP-seq uni-reads files in BED format.}
\item{chipFormat}{Location of the final allocated ChIP-seq
  output by Permseq in other selected format (tagAlign or bed).}
\item{dataNum}{The number of dataset(s) that are used.}
\item{chrList}{Chromosome list.}
\item{fragL}{Fragment length.}
\item{bowtieInfo}{Bowtie information. Includes: location of bowtie
  index, directory where the bowtie is installled, bowtie parameters values -v, -m and -p.}
\item{bwaInfo}{BWA related information. Includes: location of bwa index,
  directory where the BWA is installed, BWA parameters values for "bwa aln -n
  -o -t" and "bwa samse -n".}
\item{csemDir}{Directory where CSEM was installed.}
\item{picardDir}{Directory where PICARD jar file is saved.}
\item{outfileLoc}{Directory to store output files.}
\item{prior}{Location where the prior files are saved.}
\item{chrom.ref}{Name of the file for chromosome info. In this
file, the first and second columns are ID and size of each
chromosome, respectively.}
}
\description{
Methods included in "Prior" class: \code{names()}, \code{plot()},
\code{print()}, \code{summary()}  and \code{show()}. 
}
\examples{
showClass("Prior")
}
\author{
Ye Zheng,  M. Constanza Rojo-Alfaro.
}

