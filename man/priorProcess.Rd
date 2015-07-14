\name{priorProcess}
\alias{priorProcess}
\title{Process DNase-seq (and histone ChIP-seq) datasets}
\usage{
priorProcess(dnaseFile = NULL, histoneFile = NULL, dnaseName
= "dnase", histoneName = NULL, fragL, AllocThres
= 900, chrList = NULL, capping = 0, outfileLoc = "./",
outfile, bowtieIndex, csemDir, bowtieDir, vBowtie = 2,
mBowtie = 99, pBowtie = 8, save.files = TRUE)
}
\arguments{
\item{dnaseFile}{DNase-seq file in fastq format. For faster results, bam
  or bed files already obtained by CSEM with allocated reads can also be processed, if available.}

\item{histoneFile}{Histone ChIP-seq file in fastq format. For faster results, bam
  or bed files already obtained by CSEM with allocated reads can also be processed, if available.}

\item{dnaseName}{Name of DNase-seq data or the dataset used as DNase-seq
data in the model, default set as 'dnase'.}

\item{histoneName}{Name of histone ChIP-seq data sets. If no giving
values, histoneName is set as a vector of index number
\code{1:length(histoneFile)}.}

\item{fragL}{Fragment length.}

\item{AllocThres}{Allocation threshold. It will select reads with  scores higher than\cr
  \code{AllocThres} (allocation probability*1000). Default set at 900.}

\item{chrList}{A vector of chromosomes that will be included in
the analysis. Default set as NULL and \code{priorProcess} will get the
list from processed files. Otherwise, if given by the user, it should be consistent with the
chromosome name(s) in the corresponding fasta file(s). For more information, see details.}

\item{capping}{Maximum number of reads allowed to start at each
nucleotide position.  To avoid potential PCR amplification
artifacts, the maximum number of reads that can start at a
nucleotide position is capped.  Default is 0 (no capping).}

\item{outfileLoc}{Directory to store processed files.}

\item{outfile}{Infix of outfile name.}

\item{bowtieIndex}{Bowtie index, used in bowtie aligning}

\item{csemDir}{Directory where CSEM was installed.}

\item{bowtieDir}{Directory where Bowtie was installed}

\item{vBowtie}{Bowtie parameter. In -v mode, alignments may have
no more than vBowtie mismatches, where \code{v} may be a number from 0
through 3 set using the -v option. Default value is 2.}

\item{mBowtie}{Bowtie parameter. -m parameter instructs bowtie to
refrain from reporting any alignments for reads having more than
\code{mBowtie} reportable alignments. Default value is 99 allowing
multi-reads alignment.}

\item{pBowtie}{Bowtie parameter. The -p option causes Bowtie to launch a specified
number of parallel search threads. Each thread runs on a different
processor/core and all threads find alignments in paralle. Default value is 8.}

\item{save.files}{Option to save intermediate files created. Default set as TRUE.}

}
\value{
 A new "Prior" object is created containing the following information:

\item{dnaseName}{Name of dataset that is used as DNase-seq, especially in the "Only histone" situation, dnaseName is the selected histone ChIP-seq dataset.}
\item{dnaseAlign}{DNase-seq alignment summary information from bowtie.}
\item{dnaseKnots}{A vector of knots points for the B-spline
  functions. They are the 90,  99 and 99.9th percentiles of read
  counts.}
  \item{dnaseThres}{A vector of DNase-seq group created to generate aggregated ChIP
data. After alignment, positions which have the same DNase-seq read count are clustered into
  one group. \code{dnaseThres} is the corresponding read count number
  in each group. Each count is corresponding to one group and the grouping data 
partitions the whole genome into multiple segments.}
\item{posLoc_bychr}{Location of the files containing the group index of
  each segment of the genome.}
\item{dnaseHistone}{The input histone ChIP-seq files which are the candidates for the DNase-seq while running priorProcess procedure.}
\item{histoneName}{Name of histone ChIP-seq dataset(s). If no giving values, the histoneName would be set as a vector of index number(1:length(histoneFile)).}
\item{histoneNum}{Number of histone ChIP-seq dataset(s).}
\item{histoneAlign}{histone ChIP-seq alignment summary information from bowtie.}
\item{histoneGrpL}{The histone ChIP-seq dataset(s) selected to be included in the model after variable selection.}
\item{chipName}{Name of ChIP-seq dataset(s).}
\item{chipNum}{Number of ChIP-seq dataset(s).}
\item{chipAlign}{ChIP-seq alignment summary information from bowtie.}
\item{chipSAM}{Location of aligned ChIP-seq in SAM format.}
\item{chipAllocate}{Location of the final allocated ChIP-seq output
  by Permseq in BAM format.}
\item{chipUni}{Location of the aligned ChIP-seq uni-reads files in BED format.}
\item{chipFormat}{Location of the final allocated ChIP-seq
  output by Permseq in other selected format (tagAlign or bed).}
\item{dataNum}{Number of dataset(s) that are used.}
\item{fragL}{Fragment length.}
\item{bowtieInfo}{List of bowtie information used: bowtieIndex, bowtieDir, vBowtie, mBowtie and pBowtie.}
\item{csemDir}{Directory where CSEM was installed.}
\item{outfileLoc}{Directory to store output files.}
\item{prior}{Location where the prior files are saved.}
\item{chrom.ref}{Name of the file for chromosome info. In this
file, the first and second columns are ID and size of each
chromosome, respectively.}
}
\description{

Process DNase-seq and/or histone ChIP-seq data and construct required DNase-seq and histone ChIP-seq
files for generating priors to allocate multi-reads in ChIP-Seq data.
}
\details{
Processes DNase-seq and/or histone ChIP-seq files and generates module for further
analysis in \code{priorGenerate}. If no DNase-seq data available and do not
know which histone data could play as the DNase-seq data in the model, start from
\code{priorHistone_init} and \code{priorHistone_multi} functions.\cr

If no DNase-seq or histone ChIP-seq data available, run \code{readAllocate} directly and multi-reads will
be allocated without using prior information.\cr

If no chrList is provided,  \code{priorProcess} will generate the
list from processed files (.sam file if DNase-seq input file is in fastq
format or .bed file if DNase-seq input file is in .bam or .bed
format). Otherwise, if given by the user, it will accelerate the procedure, but the chrList should be consistent with the
chromosome name(s) in the corresponding .fa or .fasta file(s). In other
words, for example, it should be the name on the first line after ">" in .fa file.

\code{plot()}, \code{summary()}, \code{names()} and \code{print()}
methods can be used to see the information contained in "Prior"
object. To obtain the ChIP-drq alignment information from bowtie, use \code{summary()}.


}
\examples{
\dontrun{
object=priorProcess(dnaseFile=NULL, histoneFile=NULL, dnaseName='dnase',
histoneName=NULL, fragL, chrList, capping=0, outfileLoc="./",
outfile, bowtieIndex, csemDir)
}
}
\author{
Xin Zeng, M. Constanza Rojo-Alfaro, Ye Zheng
}
\keyword{methods}
\keyword{models}

