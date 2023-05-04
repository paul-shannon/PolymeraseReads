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

f.chroPlus <- system.file(package="PolymeraseReads", "extdata",
                         "chroPlus-hg38.bw")
                         # "GSM3309956_5587_5598_24205_HGC2FBGXX_J_CHR_TGACCA_R1_plus.bw")
f.chroMinus <- system.file(package="PolymeraseReads", "extdata",
                         "chroMinus-hg38.bw")
                         #  "GSM3309956_5587_5598_24205_HGC2FBGXX_J_CHR_TGACCA_R1_minus.bw")

fileList <- c(pol2=f.pol2, rna=f.rna, chroMinus=f.chroMinus, chroPlus=f.chroPlus)

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
    test_liftovers()
    test_getReads()
    test_rpkm.all.assays()
    #test_rpkm.rna()
    # test_rpkm.pol2()

} # runTests
#----------------------------------------------------------------------------------------------------
test_ctor <- function()
{
    message(sprintf("--- test_ctor"))

    geneSymbol <- "GATA2"
    stopifnot(file.exists(f.pol2))
    psr <- PolymeraseReads$new(fileList)

    checkTrue(all(c("R6", "PolymeraseReads") %in% class(psr)))
    psr$setGeneSymbol(geneSymbol)
    checkEquals(psr$getGeneSymbol(), geneSymbol)
    checkEquals(psr$getGeneID(), "2624")

    fileList <- c(pol2=f.pol2, rna=f.rna, chroPlus=f.chroPlus, chroMinus=f.chroMinus)

    psr <- PolymeraseReads$new(fileList)
    checkTrue(all(c("R6", "PolymeraseReads") %in% class(psr)))

} # test_ctor
#----------------------------------------------------------------------------------------------------
test_getTranscriptCoordinates <- function()
{
    message(sprintf("--- test_identifyTranscriptRegion"))

    geneSymbol <- "RPL13A"
    fileList <- c(pol2=f.pol2, rna=f.rna)
    psr <- PolymeraseReads$new(fileList)
    psr$setGeneSymbol(geneSymbol)

    psr$findLongestTranscript()
    txs <- psr$getTranscriptCoordinates()
    checkEquals(names(txs), c("hg38", "hg19"))
    checkEquals(as.integer(lapply(txs, start)), c(49487608, 49990865))
    checkEquals(as.integer(lapply(txs, end)), c(49493057, 49996314))
    checkEquals(as.character(seqnames(txs$hg38)), "chr19")
    checkEquals(as.character(seqnames(txs$hg19)), "chr19")
    strands <- as.character(lapply(txs, function(tx) as.character(strand(tx))))

    if(viz){
      track <- GRangesAnnotationTrack("tx", txs$hg38, color="blue")
      displayTrack(igv, track)
      }

    checkEquals(strands, c("+", "+"))


} # test_getTranscriptCoordinates
#----------------------------------------------------------------------------------------------------
test_liftovers <- function()
{
    geneSymbol <- "TAL1"
    psr <- PolymeraseReads$new(fileList)
    psr$setGeneSymbol(geneSymbol)

    psr$findLongestTranscript()
    gr.tx.hg38 <- psr$getTranscriptCoordinates()$hg38

    gr.hg19 <- psr$lift.hg38.to.hg19(gr.tx.hg38)
    gr.tx.hg38.recovered <- psr$lift.hg19.to.hg38(gr.hg19)
    checkTrue(gr.tx.hg38 == gr.tx.hg38.recovered)

} # test_liftovers
#----------------------------------------------------------------------------------------------------
test_getReads <- function()
{
    message(sprintf("--- test_getReads"))

    geneSymbol <- "TAL1"
    igv <- NA

    psr <- PolymeraseReads$new(fileList)
    psr$setGeneSymbol(geneSymbol)

    psr$findLongestTranscript(igv)
    gr.tx <- psr$getTranscriptCoordinates()$hg38

       #---------------
       # pol2 reads
       #---------------

    gr.reads.pol2 <- psr$getReads("pol2") # , igv, use.igv.roi)
    checkTrue(length(gr.reads.pol2) > 1050)

       #---------------
       # rna reads
       #---------------

    gr.reads.rna <- psr$getReads("rna", igv, use.igv.roi)
    checkTrue(length(gr.reads.rna) > 1200)
    checkTrue(length(gr.reads.rna) < 1300)

       #---------------
       # ChRO + reads
       #---------------

    gr.reads.chroPlus <- psr$getReads("chroPlus", igv, use.igv.roi)
    checkTrue(length(gr.reads.chroPlus) > 380)
    checkTrue(length(gr.reads.chroPlus) < 420)

       #---------------
       # ChRO - reads
       #---------------

    gr.reads.chroMinus <- psr$getReads("chroMinus", igv, use.igv.roi) # 1188
    checkTrue(length(gr.reads.chroMinus) > 1150)
    checkTrue(length(gr.reads.chroMinus) < 1250)

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

    psr <- PolymeraseReads$new(fileList)
    psr$setGeneSymbol(geneSymbol)

    gr.tx <- psr$getTranscriptCoordinates()
    checkTrue(width(gr.tx) > 15000)
    gr.reads <- psr$getReads("pol2")
    length(gr.reads)
    checkTrue(length(gr.reads) > 1050)
    checkTrue(length(gr.reads) < 1150)

    gr.reads <- psr$getReads("rna")
    length(gr.reads)
    checkTrue(length(gr.reads) > 1250)
    checkTrue(length(gr.reads) < 1300)

    gr.reads <- psr$getReads("chroMinus")
    length(gr.reads)
    checkTrue(length(gr.reads) > 1)
    checkTrue(length(gr.reads) < 5)

    gr.reads <- psr$getReads("chroPlus")

    tal1.rpkm <- rpkm(geneSymbol, "rna")

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
    psr <- PolymeraseReads$new(fileList)
    tal1.rpkm <- psr$rpkm(goi, assay="rna", igv=igv)
    checkEqualsNumeric(tal1.rpkm, 342, tolerance=0.2)

    goi <- "MUTYH"
    mutyh.rpkm <- psr$rpkm(goi, assay="rna")
    checkEqualsNumeric(mutyh.rpkm, 6.3, tolerance=0.2)

    goi <- "YBX1"
    goi.rpkm <- psr$rpkm(goi, assay="rna", igv=igv)
    checkEqualsNumeric(goi.rpkm, 342, tolerance=0.2)

    goi <- "RPL22"  # 3 exons, high cores
    showGenomicRegion(igv, goi)
    zoomOut(igv)
    goi.rpkm <- psr$rpkm(goi, assay="rna", igv=igv)
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

    psr <- PolymeraseReads$new(fileList)
    tal1.rpkm.0 <- psr$rpkm(goi, assay="pol2", start.site.avoidance=0, igv=igv)
    checkEqualsNumeric(tal1.rpkm.0, 6.05, tolerance=0.003)

    tal1.rpkm.500 <- psr$rpkm(goi, assay="pol2", start.site.avoidance=2000, igv=igv)
    checkEqualsNumeric(tal1.rpkm.500, 6.48, tolerance=0.003)

    goi <- "YBX1"
    showGenomicRegion(igv, goi)
    goi.rpkm.0 <- psr$rpkm(goi, assay="pol2", start.site.avoidance=0, igv=igv)
    checkEqualsNumeric(goi.rpkm.0, 1.96, tolerance=0.003)

    goi <- "YBX1"
    goi.rpkm.500 <- psr$rpkm(goi, assay="pol2", start.site.avoidance=500, igv=igv)
    checkEqualsNumeric(goi.rpkm.500, 1.49, tolerance=0.003)

}  # test_rpkm.pol2
#----------------------------------------------------------------------------------------------------
test_rpkm.all.assays <- function()
{
    message(sprintf("--- test_prkm.all.assays"))

    psr <- PolymeraseReads$new(fileList)
    goi <- "TAL1"
    psr$setGeneSymbol(goi)

        #------------------------------------------------------------
        # first, rpkm from pol2, 0, 500, and 2000 kb tss avoidance
        # note that rkpm rises slightly with each larger avoidance
        # reflecting, I think, the absence of a big tss pile-up
        # and high pol2 reads downstream
        #------------------------------------------------------------
    igv <- start.igv(goi)
    tal1.rpkm.0 <- psr$rpkm(goi, assay="pol2", start.site.avoidance=0, igv=igv)
    checkEqualsNumeric(tal1.rpkm.0, 6.05, tolerance=1e-3)

    tal1.rpkm.500 <- psr$rpkm(goi, assay="pol2", start.site.avoidance=500, igv=igv)
    checkEqualsNumeric(tal1.rpkm.500, 6.12, tolerance=1e-3)

    tal1.rpkm.2000 <- psr$rpkm(goi, assay="pol2", start.site.avoidance=2000, igv=igv)
    checkEqualsNumeric(tal1.rpkm.2000, 6.48, tolerance=1e-3)


        #-----------------
        # rna-seq
        #-----------------

    psr <- PolymeraseReads$new(fileList)
    goi <- "TAL1"
    psr$setGeneSymbol(goi)
    tal1.rna.rpkm <- psr$rpkm(goi, assay="rna", start.site.avoidance=0, igv=igv)
    checkEqualsNumeric(tal1.rna.rpkm, 34.16, tolerance=1e-3)

        #------------------------------
        # ChroPlus at for tss +0, +500
        #------------------------------

    psr <- PolymeraseReads$new(fileList)
    goi <- "TAL1"
    psr$setGeneSymbol(goi)

    tal1.chroPlus.rpkm.0 <- psr$rpkm(goi, assay="chroPlus", start.site.avoidance=0, igv=igv)
    checkEqualsNumeric(tal1.chroPlus.rpkm.0, 9.92, tolerance=1e-3)

    tal1.chroPlus.rpkm.500 <- psr$rpkm(goi, assay="chroPlus", start.site.avoidance=500, igv=igv)
    checkEqualsNumeric(tal1.chroPlus.rpkm.500, 8.95, tolerance=1e-3)

        #------------------------------
        # ChroMinus at for tss +0, +500
        #------------------------------

    psr <- PolymeraseReads$new(fileList)
    goi <- "TAL1"
    psr$setGeneSymbol(goi)

    tal1.chroMinus.rpkm.0 <- psr$rpkm(goi, assay="chroMinus", start.site.avoidance=0, igv=igv)
    checkEqualsNumeric(tal1.chroMinus.rpkm.0, 31.48, tolerance=1e-3)

    tal1.chroMinus.rpkm.500 <- psr$rpkm(goi, assay="chroMinus", start.site.avoidance=500, igv=igv)
    checkEqualsNumeric(tal1.chroMinus.rpkm.500, 30.96, tolerance=1e-3)

} # test_prkm.all.assays
#----------------------------------------------------------------------------------------------------

if(!interactive())
    runTests()
