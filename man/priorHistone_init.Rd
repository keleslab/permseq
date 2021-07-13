\name{priorHistone_init}
\alias{priorHistone_init}
\title{Process histone ChIP-seq dataset to build the prior without DNase-seq}
\usage{
priorHistone_init(histoneFile = NULL, histoneName = NULL, chipFile = NULL, fragL = 200,
AllocThres = 900, chrList, capping = 0, outfileLoc = "./", 
bowtieDir, bowtieIndex, vBowtie = 2, mBowtie = 99,
pBowtie = 8, bwaDir, bwaIndex, nBWA = 2, oBWA = 1, tBWA = 8, mBWA = 99,
csemDir, picardDir, saveFiles = TRUE)
}
\arguments{
\item{histoneFile}{Input histone ChIP-seq file in fastq format. For faster results, sam
  formatted file after alignment including multi-mapping reads or bam
  or bed files already obtained by CSEM with allocated reads can also be
  processed, if available. Otherwise, it is better start from the fastq
  formatted file. Default value is NULL.}

\item{histoneName}{Name of the histone ChIP-seq data. If multiple Histone
data are used, it needs to be formed into a vector.}

\item{chipFile}{Input ChIP-seq files, in fastq format or sam format to
  save time if it is already aligned and includes multi-mapping
  reads. The default value will be NULL.}

\item{fragL}{Average fragment length. The default value is 200.}

\item{AllocThres}{Allocation threshold. It will select reads with scores higher than \code{AllocThres} (allocation probability*1000). 
Default set at 900.}

\item{chrList}{A vector of chromosomes that will be included in
the analysis. Default set as NULL and \code{priorHistone_init} will get the
list from processed files. Otherwise, if given by the user, it should be consistent with the
chromosome name(s) in the corresponding fasta file(s). For more information, see details.}

\item{capping}{Maximum number of reads allowed to start at each
nucleotide position.  To avoid potential PCR amplification
artifacts, the maximum number of reads that can start at a
nucleotide position is capped.  Default is 0 (no capping, i.e. no maximum restriction).}

\item{outfileLoc}{Directory to store processed files. Default is set to "./".}


\item{bowtieDir}{Directory where Bowtie was installed. Default will be NULL.}

\item{bowtieIndex}{Bowtie index, used in bowtie alignment. Users can
  select the aligner, Bowtie or BWA, by providing the corresponding
  index. Default will be NULL.}

\item{vBowtie}{Bowtie parameter. In -v mode, alignments may have
no more than vBowtie mismatches, where \code{v} may be a number from 0
through 3 set using the -v option. Default value is 2.}

\item{mBowtie}{Bowtie parameter. -m parameter instructs bowtie to
refrain from reporting any alignments for reads having more than
\code{mBowtie} reportable alignments. Default value is 99 allowing
all kinds of multi-reads alignment.}

\item{pBowtie}{Bowtie parameter. The -p option causes Bowtie to launch a specified
number of parallel search threads. Each thread runs on a different
processor/core and all threads find alignments in paralle. Default value is 8.}

\item{bwaDir}{Directory where BWA was installed. Default set as NULL.}

\item{bwaIndex}{BWA index used in BWA alignment. Users can specify the aligner, Bowtie or
  BWA, by specifying the index that will be used. Default set as NULL.}

\item{nBWA}{BWA paramter. In "bwa aln -n" mode, if it is an integer, it denotes the
  maximum edit distances including mismatch and gap open. Otherwise, it
  will be the fraction of missing alignments given 2\% uniform base errr
  rate. Default value is 2.}

\item{oBWA}{BWA parameter. In "bwa aln -o" mode, it specifies the
  maximum number of gap open. Default set as 1.}

\item{tBWA}{BWA parameter. In "bwa aln -t" mode, it is the number of threads in
  multi-threading mode. Default set as 8.}

\item{mBWA}{BWA parameter. In "bwa samse -n", it restricts the maximum number of
  alignments to output for each read. If a read has more hits, the XA
  tag will not be written. Default set as 99.}


\item{csemDir}{Directory where CSEM was installed. The default value is NULL.}
\item{picardDir}{Directory where PICARD jar file is saved.}

\item{saveFiles}{Option to save intermediate files created. Default set as TRUE.}

}

\value{
  A new "Prior" object is created.

  First, for each histone ChIP-seq dataset within \code{dnaseHistone}

  \item{dnaseKnots}{Knots for B-spline functions. They are the 90, 99
  and 99.9th percentiles of read counts.}

  \item{dnaseThres}{A vector of DNase-seq group created to generate aggregated ChIP
data. After alignment, positions which have the same DNase-seq read count are clustered into
  one group. \code{dnaseThres} is the corresponding read count number
  in each group. Each count is corresponding to one group and the grouping data 
partitions the whole genome into multiple segments.}

\item{posLoc_bychr}{Location of the files containing the group index of
  each segment of the genome.}


\item{posLoc3_bychr}{Location of the files containing which segments of
  the genome are in which group based on the trinary Histone positions
  according to 90 and 99th percentiles of read counts.}


Other elements of the "Prior" include:

\item{chipName}{Name of ChIP-seq dataset(s).}
\item{chipNum}{Number of ChIP-seq dataset(s).}
\item{chipAlign}{ChIP-seq alignment summary information from bowtie.}
\item{chipSAM}{Location of aligned ChIP-seq in SAM format.}
\item{chipUni}{Location of the aligned ChIP-seq uni-reads files in BED format.}
\item{histoneName}{Name of histone ChIP-seq dataset(s). If no giving values, histoneName is set as a vector of index number(1:length(histoneFile)).}
\item{histoneNum}{Number of Histone dataset(s).}
\item{histoneAlign}{Histone alignment summary information from bowtie.}
\item{chrList}{Chromosome list.}
\item{fragL}{Fragment length, given as a parameter.}
\item{bowtieInfo}{List of bowtie information used: bowtieIndex,
  bowtieDir, vBowtie, mBowtie and pBowtie.}
\item{bwaInfo}{List of BWA related information: bwaDir, bwaIndex, nBWA,
  oBWA, tBWA, mBWA.}
\item{csemDir}{Directory of CSEM given as a parameter.}
\item{picardDir}{Directory where PICARD jar file is saved.}
\item{chrom.ref}{Name of the file for chromosome info, given as a parameter.}
\item{outfileLoc}{Directory where processed files are (given as an argument).}


}
\description{
If no DNase-seq dataset available, process Histone data and select the
histone ChIP-seq dataset who shows the strongest relationship with
ChIP-seq data. The selected histone ChIP-seq data will be used as DNase-seq
data in the model.
}
\details{
Process histone ChIP-seq files and generate module for further analysis in \cr
\code{priorHistone_multi}.\cr

This function processes the histone ChIP-seq files and generates a marginal plot
(marginal_plot.pdf) stored in \code{outfileLoc} to decide which histone
data should be used as DNase-seq data (choose the histone ChIP-seq data which has the
most increasing relationship with ChIP data). \cr

If no chrList is provided,  \code{priorHistone_init} will generate the
list from processed files (.sam file if Histone input file is in fastq
format or .bed file if Histone input file is in .bam or .bed
format). Otherwise, if given by the user, it will accelerate the procedure, but the chrList should be consistent with the chromosome name(s)
in the corresponding .fa or .fasta file(s). In other words, for example, it should be the name on the first line after ">" in .fa file.

Users can select from Bowtie and BWA to do the alignment by providing
the corresponding index and leaving the other as default value NULL. If
both indices are provided, the package will automatically use Bowtie to do the multi-mapping reads alignment. 

Aligned sam file will go through filtering process to remove
duplicates. By default, 'samtools rmdup -s' will be used to carry out
such function. If user provides the PICARD jar path, PICARD will be
used.

\code{plot()}, \code{summary()}, \code{names()} and \code{print()}
methods can be used to see the information contained in "Prior"
object. To obtain the ChIP-seq (or Histone) alignment information from bowtie, use
\code{summary()}.

}
\examples{
\dontrun{
object = priorHistone_init(histoneFile = NULL, histoneName = NULL,
  chipFile, fragL, AllocThres = 900, chrList, capping = 0, outfileLoc = "./",
  bowtieDir,  bowtieIndex,  vBowtie = 2,
  mBowtie = 99, pBowtie = 8, csemDir, picardDir, saveFiles = TRUE)
}
}
\author{
Xin Zeng,  M. Constanza Rojo-Alfaro, Ye Zheng
}
\keyword{methods}
\keyword{models}

