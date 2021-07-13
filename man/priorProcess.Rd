\name{priorProcess}
\alias{priorProcess}
\title{Process DNase-seq (and histone ChIP-seq) datasets}
\usage{
priorProcess(dnaseFile = NULL, histoneFile = NULL, dnaseName
= "dnase", histoneName = NULL, fragL = 200, AllocThres
= 900, chrList = NULL, capping = 0, outfileLoc = "./",
outfile = "dnase_histone", bowtieDir, bowtieIndex, vBowtie = 2,
mBowtie = 99, pBowtie = 8, bwaDir, bwaIndex, nBWA = 2, oBWA = 1, tBWA =
8, mBWA = 99, csemDir, picardDir, chrom.ref, saveFiles = TRUE)
}
\arguments{
\item{dnaseFile}{DNase-seq file in fastq format. For faster results, sam
  formatted file after alignment including multi-mapping reads or bam
  or bed files already obtained by CSEM with allocated reads can also be
  processed, if available. Otherwise, it is better start from the fastq
  formatted file. The default value is NULL.}

\item{histoneFile}{Histone ChIP-seq file in fastq format. For faster results, sam
  formatted file after alignment including multi-mapping reads or bam
  or bed files already obtained by CSEM with allocated reads can also be
  processed, if available. Otherwise, it is better start from the fastq
  formatted file. Default value is set to NULL.}

\item{dnaseName}{Name of DNase-seq data or the dataset used as DNase-seq
data in the model, default set as 'dnase'.}

\item{histoneName}{Name of histone ChIP-seq data sets. If no giving
values, histoneName is set as a vector of index number
\code{1:length(histoneFile)}.}

\item{fragL}{Average fragment length with default value 200.}

\item{AllocThres}{Allocation threshold. It will select reads with  scores higher than\cr
  \code{AllocThres} (allocation probability*1000). Default set at 900.}

\item{chrList}{A vector of chromosomes that will be included in
the analysis. Default set as NULL and \code{priorProcess} will get the
list from processed files. Otherwise, if given by the user, it should be consistent with the
chromosome name(s) in the corresponding fasta file(s). For more information, see details.}
 
\item{capping}{ Maximum number of reads allowed at each
nucleotide position.  To avoid potential PCR amplification
artifacts, the maximum number of reads that can start at a
nucleotide position is capped.  Default is 0 (no capping, i.e. no maximum restriction).}

\item{outfileLoc}{Directory to store processed files.}

\item{outfile}{Infix of outfile name. Default set as "dnase_histone"
  indicating the prior is constructed using DNase-seq and Histone data.}


\item{bowtieDir}{Directory where Bowtie was installed, default set as NULL.}


\item{bowtieIndex}{Bowtie index, used in bowtie aligning. Default value
  is NULL and users can specify the selection of aligner, Bowtie or BWA,
  by providing the corresponding index.}

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

\item{csemDir}{Directory where CSEM was installed.}

\item{picardDir}{Directory where PICARD jar file is saved.}

\item{chrom.ref}{reference genome index summary information. First line
  is the number of chromosomes in the index, either bwaIndex or
  bowtieIndex, including chrM. Second line is the size of each
  chromosome. Third line is the name of chromsome.}

\item{saveFiles}{Option to save intermediate files created. Default set as TRUE.}

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
\item{histoneName}{Name of histone ChIP-seq dataset(s). If no giving values, the histoneName would be set as a vector of index number(1:length(histoneFile)).}
\item{histoneNum}{Number of histone ChIP-seq dataset(s).}
\item{histoneAlign}{histone ChIP-seq alignment summary information from bowtie.}
\item{dataNum}{Number of dataset(s) that are used.}
\item{chrList}{Chromosome list.}
\item{fragL}{Fragment length.}
\item{bowtieInfo}{List of Bowtie related information: bowtieDir,
  bowtieIndex, vBowtie, mBowtie and pBowtie.}
\item{bwaInfo}{List of BWA related information: bwaDir, bwaIndex, nBWA,
  oBWA, tBWA, mBWA.}
\item{csemDir}{Directory where CSEM was installed.}
\item{picardDir}{Directory where PICARD jar file is saved.}
\item{outfileLoc}{Directory to store output files.}
\item{chrom.ref}{Name of the file for chromosome info.}
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

Users can select from Bowtie and BWA to do the alignment by providing
the corresponding index and leaving the other as default value NULL. If
both indices are provided, the package will automatically use Bowtie to do the multi-mapping reads alignment. 

DNase-seq and/or Histone aligned sam file will go through filtering
process to remove duplicates. By default, 'samtools rmdup -s' function
will be used. PICARD jar can take over if PICARD jar path is provided.

\code{plot()}, \code{summary()}, \code{names()} and \code{print()}
methods can be used to see the information contained in "Prior"
object. To obtain the ChIP-drq alignment information from bowtie, use \code{summary()}.


}
\examples{
\dontrun{
object = priorProcess(dnaseFile = NULL, histoneFile = NULL, dnaseName = 'dnase',
histoneName = NULL, fragL = 200, chrList, capping = 0, outfileLoc = "./",
outfile, bowtieIndex, csemDir, picardDir, chrom.ref, saveFiles = FALSE)
}
}
\author{
Xin Zeng, M. Constanza Rojo-Alfaro, Ye Zheng
}
\keyword{methods}
\keyword{models}

