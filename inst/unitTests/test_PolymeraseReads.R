library(RUnit)
library(PolymeraseReads)
library(GenomicRanges)
library(rtracklayer)
#----------------------------------------------------------------------------------------------------
chr19.de.genes <- c("RPL131", "RPL18", "RPL18A", "EEF2", "PTB1", "BSG", "RPS28")
goi <- chr19.de.genes[1]
viz <- FALSE
if(exists("viz") && viz){
    if(!exists("igv")){
        igv <<- start.igv(goi, "hg38")
    } else {
        showGenomicRegion(igv, goi)
        }
    } # viz
#----------------------------------------------------------------------------------------------------
f.pol2 <- system.file(package="PolymeraseReads", "extdata",
                      "NS.1828.002.N704---N506.MB_OHRI_Jurkat_POL2_CI_CPM.bw")
f.rna <- system.file(package="PolymeraseReads", "extdata",
                     "RNAseq_IRCM-1793_Hs_Jurkat_DoxN_Rep1_CPM.bw")
#----------------------------------------------------------------------------------------------------
display.rna.chr1 <- function()
{
    igv <- start.igv("TAL1", "hg38")
    gr.roi <- GRanges(seqnames="chr1", IRanges(start=1, end=250000000))
    gr <- import(f.rna, which=gr.roi, format="bigwig")
    track <- GRangesQuantitativeTrack("RNA", gr, color="random", autoscale=TRUE)
    displayTrack(igv, track)

} # display.rna.chr1
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
    test_ctor()
    test_getTranscriptCoordinates()
    test_getReads()
    test_rpkm.rna()
    test_rpkm.pol2()

} # runTests
#----------------------------------------------------------------------------------------------------
test_ctor <- function()
{
    message(sprintf("--- test_ctor"))

    geneSymbol <- "GATA2"
    stopifnot(file.exists(f.pol2))
    psr <- PolymeraseReads$new(geneSymbol, f.pol2)
    checkTrue(all(c("R6", "PolymeraseReads") %in% class(psr)))
    checkEquals(psr$getGeneSymbol(), geneSymbol)
    checkEquals(psr$getGeneID(), "2624")

    psr <- PolymeraseReads$new(geneSymbol, f.pol2, f.rna)
    checkTrue(all(c("R6", "PolymeraseReads") %in% class(psr)))
    checkEquals(psr$getGeneSymbol(), geneSymbol)
    checkEquals(psr$getGeneID(), "2624")

} # test_ctor
#----------------------------------------------------------------------------------------------------
test_getTranscriptCoordinates <- function()
{
    message(sprintf("--- test_identifyTranscriptRegion"))

    geneSymbol <- "RPL13A"
    psr <- PolymeraseReads$new(geneSymbol, f.pol2, f.rna)

    gr.tx <- psr$getTranscriptCoordinates()
    tbl.tx <- as.data.frame(gr.tx)
    with(tbl.tx, checkEquals(as.character(seqnames), "chr19"))
    checkEqualsNumeric(width(gr.tx), 5450, tol=1000)
    with(tbl.tx, checkEquals(as.character(strand), "+"))

    if(viz){
      track <- GRangesAnnotationTrack("tx", gr.tx, color="blue")
      displayTrack(igv, track)
      }

} # test_getTranscriptCoordinates
#----------------------------------------------------------------------------------------------------
test_getReads <- function()
{
    message(sprintf("--- test_getReads"))

    geneSymbol <- "RPL13A"
    psr <- PolymeraseReads$new(geneSymbol, f.pol2, f.rna)

    gr.tx <- psr$getTranscriptCoordinates()
    gr.reads <- psr$getReads()
    checkTrue(length(gr.reads) > 400)
    checkTrue(length(gr.reads) < 500)

    if(viz){
       track <- GRangesQuantitativeTrack("pol2 reads", gr.reads, color="darkgreen")
       displayTrack(igv, track)
       roi <- getGenomicRegion(igv)
       gr.roi <- GRanges(seqnames=roi$chrom, IRanges(start=roi$start, end=roi$end))
       gr.bw <- import(f.rna, which=gr.roi, format="bigwig")
         # most of the reads have very low scores, presumed to be noise
       hist(gr.bw$score)
       gr.bw.trimmed <- gr.bw[gr.bw$score > 1]
       hist(gr.bw.trimmed$score)
       total.reads <- round(sum(gr.bw.trimmed$score)) # 442441 on chr3
       track <- GRangesQuantitativeTrack("rna-1", gr.bw.trimmed, autoscale=TRUE, min=0.9, max=300, color="brown")
       displayTrack(igv, track)
       }

} # test_getReads
#----------------------------------------------------------------------------------------------------
# show rna-seq and pol2/ser2 across the entire chromosome, to manually
# identify expressed genes
display.chromosome <- function()
{
   gr.roi <- GRanges(seqnames="chr19", IRanges(1,58617616))

   gr.bw <- import(f.pol2, which=gr.roi, format="bigwig")
      # most of the reads have very low scores, presumed to be noise
   fivenum(gr.bw$score)
   gr.bw.trimmed <- gr.bw[gr.bw$score > 1]
   total.reads <- round(sum(gr.bw.trimmed$score)) # 442441 on chr3
   title <- "chr19 pol2"
   track <- GRangesQuantitativeTrack(title,gr.bw.trimmed, autoscale=TRUE)
   displayTrack(igv, track)


      # most of the reads have very low scores, presumed to be noise
   gr.bw <- import(f.rna, which=gr.roi, format="bigwig")
   hist(gr.bw$score)
   gr.bw.trimmed <- gr.bw[gr.bw$score > 1]
   hist(gr.bw.trimmed$score)
   total.reads <- round(sum(gr.bw.trimmed$score)) # 442441 on chr3
   track <- GRangesQuantitativeTrack("chr19 rna", gr.bw.trimmed, autoscale=TRUE, min=0.9, max=300, color="brown")
   displayTrack(igv, track)

} # display.chromosome
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

   roi <- getGenomicRegion(igv)
   gr.roi <- with(roi, GRanges(seqnames=chrom, IRanges(start, end)))
   gr.bw <- import(f.pol2, which=gr.roi, format="bigwig")
   # most of the reads have very low scores, presumed to be noise
   fivenum(gr.bw$score)
   gr.bw.trimmed <- gr.bw[gr.bw$score > 1]
   total.reads <- round(sum(gr.bw.trimmed$score)) # 442441 on chr3
   title <- sprintf("%s pol2", gene)
   track <- GRangesQuantitativeTrack(title,gr.bw.trimmed, autoscale=TRUE)
   displayTrack(igv, track)

} # display.gene
#----------------------------------------------------------------------------------------------------
test_tal1 <- function()
{
    message(sprintf("--- test_tal1"))

    geneSymbol <- "TAL1"
    psr <- PolymeraseReads$new(geneSymbol, f.pol2, f.rna)

    gr.tx <- psr$getTranscriptCoordinates()
    gr.reads <- psr$getReads()
    checkTrue(length(gr.reads) > 1050)
    checkTrue(length(gr.reads) < 1150)
    tal1.rpkm <- rpkm(geneSymbol, f.rna)

    if(viz){
       track <- GRangesQuantitativeTrack("pol2 reads", gr.reads, color="darkgreen")
       displayTrack(igv, track)
       roi <- getGenomicRegion(igv)
       gr.roi <- GRanges(seqnames=roi$chrom, IRanges(start=roi$start, end=roi$end))
       gr.bw <- import(f.rna, which=gr.roi, format="bigwig")
         # most of the reads have very low scores, presumed to be noise
       hist(gr.bw$score)
       gr.bw.trimmed <- gr.bw[gr.bw$score > 1]
       hist(gr.bw.trimmed$score)
       total.reads <- round(sum(gr.bw.trimmed$score)) # 442441 on chr3
       track <- GRangesQuantitativeTrack("rna-1", gr.bw.trimmed, autoscale=TRUE, min=0.9, max=300, color="brown")
       displayTrack(igv, track)
       }

} # test_tal1
#----------------------------------------------------------------------------------------------------
test_rpkm.rna <- function()
{
    message(sprintf("--- test_rpkm.rna"))

    goi <- "TAL1"
    psr <- PolymeraseReads$new(goi, f.pol2, f.rna)
    tal1.rpkm <- psr$rpkm(goi, assay="rna", igv=igv)
    checkEqualsNumeric(tal1.rpkm, 34, tolerance=0.2)


    goi <- "MUTYH"
    mutyh.rpkm <- psr$rpkm(goi, assay="rna")
    checkEqualsNumeric(mutyh.rpkm, 6.3, tolerance=0.2)

    goi <- "YBX1"
    goi.rpkm <- psr$rpkm(goi, assay="rna", igv=igv)
    checkEqualsNumeric(goi.rpkm, 342, tolerance=0.2)

    goi <- "RPL22"  # 3 exons, high cores
    goi.rpkm <- psr$rpkm(goi, assay="rna")
    checkEqualsNumeric(goi.rpkm, 264, tolerance=0.2)

} # test_rpkm.rna
#----------------------------------------------------------------------------------------------------
test_rpkm.pol2 <- function()
{
    message(sprintf("--- test_rpkm.pol2"))

    goi <- "TAL1"
    if(FALSE){
        igv <- start.igv(goi, "hg38")
        }

    psr <- PolymeraseReads$new(goi, f.pol2, f.rna)
    tal1.rpkm.0 <- psr$rpkm(goi, assay="pol2", start.site.avoidance=0, igv=igv)
    checkEqualsNumeric(tal1.rpkm.0, 6.05, tolerance=0.003)

    tal1.rpkm.500 <- psr$rpkm(goi, assay="pol2", start.site.avoidance=2000, igv=igv)
    checkEqualsNumeric(tal1.rpkm.500, 6.122, tolerance=0.003)

    goi <- "YBX1"
    showGenomicRegion(igv, goi)
    goi.rpkm.0 <- psr$rpkm(goi, assay="pol2", start.site.avoidance=0, igv=igv)
    checkEqualsNumeric(goi.rpkm.0, 1.96, tolerance=0.003)

    goi <- "YBX1"
    goi.rpkm.500 <- psr$rpkm(goi, assay="pol2", start.site.avoidance=500, igv=igv)
    checkEqualsNumeric(goi.rpkm.500, 1.49, tolerance=0.003)

}  # test_rpkm.pol2
#----------------------------------------------------------------------------------------------------

if(!interactive())
    runTests()
