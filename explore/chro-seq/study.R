library(RUnit)
library(PolymeraseReads)
library(GenomicRanges)
library(rtracklayer)
#----------------------------------------------------------------------------------------------------
goi <- "MIR181A1"
if(!exists("igv")){
    igv <- start.igv(goi, "hg19")
    zoomOut(igv)
    }
#----------------------------------------------------------------------------------------------------
f.pol2 <- system.file(package="PolymeraseReads", "extdata",
                      "NS.1828.002.N704---N506.MB_OHRI_Jurkat_POL2_CI_CPM.bw")
f.rna <- system.file(package="PolymeraseReads", "extdata",
                     "RNAseq_IRCM-1793_Hs_Jurkat_DoxN_Rep1_CPM.bw")
#----------------------------------------------------------------------------------------------------
# show the longest transcript
display.tx <- function(gene)
{
   psr <- PolymeraseReads$new(gene, f.pol2, f.rna)
   gr.tx <- psr$getTranscriptCoordinates()
   tbl.tx <- as.data.frame(gr.tx)
   title <- sprintf("%s tx", gene)
   track <- GRangesAnnotationTrack(title, gr.tx, color="darkgreen")
   displayTrack(igv, track)

} # display.tx
#----------------------------------------------------------------------------------------------------
pol.2 <- function()
{
   roi <- getGenomicRegion(igv)
   gr.roi <- with(roi, GRanges(seqnames=chrom, IRanges(start, end)))

   # gr.roi <- GRanges(seqnames="chr19", IRanges(1, 58617616))
   gr.bw <- import(f.pol2, which=gr.roi, format="bigwig")
   # most of the reads have very low scores, presumed to be noise
   fivenum(gr.bw$score)
   gr.bw.trimmed <- gr.bw[gr.bw$score > 0]
   total.reads <- round(sum(gr.bw.trimmed$score)) # 442441 on chr3
   title <- "pol2"
   track <- GRangesQuantitativeTrack(title,gr.bw.trimmed, autoscale=FALSE, min=0, max=42)
   displayTrack(igv, track)

} # pol.2
#----------------------------------------------------------------------------------------------------
# show rna-seq across the entire chromosome, to manually identify expressed genes
display.gene.pol2 <- function(gene)
{
   showGenomicRegion(igv, gene)
   zoomOut(igv)

   removeTracksByName(igv,  getTrackNames(igv)[-(1:2)])

   psr <- PolymeraseReads$new(gene, f.pol2, f.rna)
   gr.tx <- psr$getTranscriptCoordinates()
   tbl.tx <- as.data.frame(gr.tx)
   title <- sprintf("%s tx", gene)
   track <- GRangesAnnotationTrack(title, gr.tx, color="darkgreen")
   displayTrack(igv, track)


} # display.gene
#----------------------------------------------------------------------------------------------------
rna.seq <- function()
{
   roi <- getGenomicRegion(igv)
   gr.roi <- GRanges(seqnames="chr19", IRanges(1, 58617616))
   gr.bw <- import(f.rna, which=gr.roi, format="bigwig")
         # most of the reads have very low scores, presumed to be noise
   #hist(gr.bw$score)
   #gr.bw.trimmed <- gr.bw[gr.bw$score > 1]
   #    hist(gr.bw.trimmed$score)
   total.reads <- round(sum(gr.bw$score)) # 442441 on chr3
   track <- GRangesQuantitativeTrack("rna-1", gr.bw, autoscale=FALSE, min=0.9, max=300, color="brown")
   displayTrack(igv, track)

} # rna.seq
#----------------------------------------------------------------------------------------------------
# aligned to hg19
chro.seq <- function()
{
   roi <- getGenomicRegion(igv)
   gr.roi <- with(roi, GRanges(seqnames=chrom, IRanges(start, end)))
   gr.roi <- GRanges(seqnames="chr1", IRanges(1, 249250621))
   data.dir <- "~/github/PolymeraseReads/inst/extdata"
   f.minus <- "GSM3309956_5587_5598_24205_HGC2FBGXX_J_CHR_TGACCA_R1_minus.bw"
   f.plus  <- "GSM3309956_5587_5598_24205_HGC2FBGXX_J_CHR_TGACCA_R1_plus.bw"

   gr.minus.bw <- import(file.path(data.dir, f.minus), which=gr.roi, format="bigwig")
   #gr.minus.bw$score <- -1 * gr.minus.bw$score
   gr.plus.bw <- import(file.path(data.dir, f.plus), which=gr.roi, format="bigwig")

   #export(gr.minus.bw, con="farsa-minus.bw", format="bigwig")
   #export(gr.plus.bw, con="farsa-plus.bw", format="bigwig")

   title <- "ChRO +"
   limit <- 500
   track <- GRangesQuantitativeTrack(title, gr.plus.bw, autoscale=FALSE, min=0, max=limit,
                                     color="red")
   displayTrack(igv, track)

   title <- "ChRO -"
   track <- GRangesQuantitativeTrack(title, gr.minus.bw, autoscale=FALSE, min=(-1*limit), max=0,
                                     color="blue")
   displayTrack(igv, track)

   showGenomicRegion(igv, goi)
   showGenomicRegion(igv, "chr1:1-249250621")

}  # chro.seq
#----------------------------------------------------------------------------------------------------
