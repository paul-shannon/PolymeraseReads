% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PolymeraseReads.R
\name{PolymeraseReads}
\alias{PolymeraseReads}
\title{PolymeraseReads}
\description{
calculates rpkm of pol2-ser2 binding from bigwig data, omitting first exons
}
\examples{
  psr <- PolymeraseReads$new(geneSymbols="GATA2")

}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-PolymeraseReads-new}{\code{PolymeraseReads$new()}}
\item \href{#method-PolymeraseReads-setGeneSymbol}{\code{PolymeraseReads$setGeneSymbol()}}
\item \href{#method-PolymeraseReads-getGeneSymbol}{\code{PolymeraseReads$getGeneSymbol()}}
\item \href{#method-PolymeraseReads-getGeneID}{\code{PolymeraseReads$getGeneID()}}
\item \href{#method-PolymeraseReads-getStrand}{\code{PolymeraseReads$getStrand()}}
\item \href{#method-PolymeraseReads-findLongestTranscript}{\code{PolymeraseReads$findLongestTranscript()}}
\item \href{#method-PolymeraseReads-getTranscriptCoordinates}{\code{PolymeraseReads$getTranscriptCoordinates()}}
\item \href{#method-PolymeraseReads-lift.hg19.to.hg38}{\code{PolymeraseReads$lift.hg19.to.hg38()}}
\item \href{#method-PolymeraseReads-lift.hg38.to.hg19}{\code{PolymeraseReads$lift.hg38.to.hg19()}}
\item \href{#method-PolymeraseReads-getReads}{\code{PolymeraseReads$getReads()}}
\item \href{#method-PolymeraseReads-rpkm}{\code{PolymeraseReads$rpkm()}}
\item \href{#method-PolymeraseReads-clone}{\code{PolymeraseReads$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PolymeraseReads-new"></a>}}
\if{latex}{\out{\hypertarget{method-PolymeraseReads-new}{}}}
\subsection{Method \code{new()}}{
Creates a new instance of this [R6][R6::R6Class] class.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolymeraseReads$new(fileList, txdb)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{fileList}}{named character vector, a list of bigwig files where
assays (names) can be one or more of rna, pol2, chroPlus, chroMinus}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
a new instance of PolymeraseReads
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PolymeraseReads-setGeneSymbol"></a>}}
\if{latex}{\out{\hypertarget{method-PolymeraseReads-setGeneSymbol}{}}}
\subsection{Method \code{setGeneSymbol()}}{
accessor for the object's geneSybmol field
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolymeraseReads$setGeneSymbol(geneSymbol)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{geneSymbol}}{character, the HUGO symbol of the gene if interest}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PolymeraseReads-getGeneSymbol"></a>}}
\if{latex}{\out{\hypertarget{method-PolymeraseReads-getGeneSymbol}{}}}
\subsection{Method \code{getGeneSymbol()}}{
accessor for the object's geneSybmol field
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolymeraseReads$getGeneSymbol()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
the geneSymbol
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PolymeraseReads-getGeneID"></a>}}
\if{latex}{\out{\hypertarget{method-PolymeraseReads-getGeneID}{}}}
\subsection{Method \code{getGeneID()}}{
accessor for the object's entrez geneID
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolymeraseReads$getGeneID()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
the geneID
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PolymeraseReads-getStrand"></a>}}
\if{latex}{\out{\hypertarget{method-PolymeraseReads-getStrand}{}}}
\subsection{Method \code{getStrand()}}{
accessor for the strand of the longest tx
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolymeraseReads$getStrand()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
character, + or - (do not expect "*")
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PolymeraseReads-findLongestTranscript"></a>}}
\if{latex}{\out{\hypertarget{method-PolymeraseReads-findLongestTranscript}{}}}
\subsection{Method \code{findLongestTranscript()}}{
find the coordinates of the longest transcript
  in both hg38 and hg19
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolymeraseReads$findLongestTranscript(igv = NA)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{igv}}{igvR instance, default NA}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
a GRanges object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PolymeraseReads-getTranscriptCoordinates"></a>}}
\if{latex}{\out{\hypertarget{method-PolymeraseReads-getTranscriptCoordinates}{}}}
\subsection{Method \code{getTranscriptCoordinates()}}{
return the coordinates of the longest transcript
  in both hg38 and hg19
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolymeraseReads$getTranscriptCoordinates()}\if{html}{\out{</div>}}
}

\subsection{Returns}{
a list of 2 GRanges objects
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PolymeraseReads-lift.hg19.to.hg38"></a>}}
\if{latex}{\out{\hypertarget{method-PolymeraseReads-lift.hg19.to.hg38}{}}}
\subsection{Method \code{lift.hg19.to.hg38()}}{
liftover GRanges from hg19 to hg38
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolymeraseReads$lift.hg19.to.hg38(gr)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{gr}}{GRanges, the hg19 data structure}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
GRanges
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PolymeraseReads-lift.hg38.to.hg19"></a>}}
\if{latex}{\out{\hypertarget{method-PolymeraseReads-lift.hg38.to.hg19}{}}}
\subsection{Method \code{lift.hg38.to.hg19()}}{
liftover GRanges from hg38 to hg18
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolymeraseReads$lift.hg38.to.hg19(gr)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{gr}}{GRanges, the hg38 data structure}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
GRanges
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PolymeraseReads-getReads"></a>}}
\if{latex}{\out{\hypertarget{method-PolymeraseReads-getReads}{}}}
\subsection{Method \code{getReads()}}{
the GRanges of scored and aligned reads of
the longest transcript of geneSymbol
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolymeraseReads$getReads(assay, igv = NA, use.igv.roi = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{assay}}{character, the name of a single assay}

\item{\code{igv}}{igvR instance, default NA}

\item{\code{use.igv.roi}}{logical query igv for current region of interest}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
GRanges
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PolymeraseReads-rpkm"></a>}}
\if{latex}{\out{\hypertarget{method-PolymeraseReads-rpkm}{}}}
\subsection{Method \code{rpkm()}}{
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolymeraseReads$rpkm(
  geneSymbol,
  score.threshold = 0,
  assay,
  start.site.avoidance = 500,
  igv = NA
)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{geneSymbol}}{character HUGO gene symbol name}

\item{\code{score.threshold}}{numeric ignore reads scoring less than this value}

\item{\code{assay}}{character either "rna" or "pol2"}

\item{\code{start.site.avoidance}}{numeric default 500bp, for pol2 only, permits
rpkm calculation to avoid distortion by}

\item{\code{igv}}{default NA, a reference to an igvR instance}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
the rpkm score
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-PolymeraseReads-clone"></a>}}
\if{latex}{\out{\hypertarget{method-PolymeraseReads-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{PolymeraseReads$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
