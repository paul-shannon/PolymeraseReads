library(biomaRt)
library(PolymeraseReads)

dataset <- "hsapiens_gene_ensembl"
hg38.mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset=dataset)

attribs <- c("chromosome_name", "transcription_start_site", "hgnc_symbol", "strand")
if(!exists("tbl.chr1")){
   tbl.chr1 <- getBM(attributes=attribs,
                     filters="chromosome_name",
                     value="1",
                     mart=hg38.mart)
   dim(tbl.chr1)
   geneSymbols <- sort(unique(tbl.chr1$hgnc_symbol))
   deleters <- which(nchar(geneSymbols) == 0)
   length(deleters)
   geneSymbols <- geneSymbols[-deleters]
   length(geneSymbols)
   head(geneSymbols)
   }
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
if(!exists("psr"))
    psr <- PolymeraseReads$new(fileList)

#----------------------------------------------------------------------------------------------------
rpkm.4way <- function(goi)
{
  timeStamp <- format(Sys.time(), "%H:%M:%S")

  psr$setGeneSymbol(goi)
  tx <- psr$findLongestTranscript()

  rna.rpkm           <- psr$rpkm(goi, assay="rna", start.site.avoidance=0)
  pol2.rpkm.0        <- psr$rpkm(goi, assay="pol2", start.site.avoidance=0)
  pol2.rpkm.500      <- psr$rpkm(goi, assay="pol2", start.site.avoidance=500)
  chroPlus.rpkm.0    <- psr$rpkm(goi, assay="chroPlus", start.site.avoidance=0)
  chroPlus.rpkm.500  <- psr$rpkm(goi, assay="chroPlus", start.site.avoidance=500)
  chroMinus.rpkm.0   <- psr$rpkm(goi, assay="chroMinus", start.site.avoidance=0)
  chroMinus.rpkm.500 <- psr$rpkm(goi, assay="chroMinus", start.site.avoidance=500)

  message(sprintf("--- rpkm.4way %s: %s %5.2f", goi, timeStamp, rna.rpkm))

  data.frame(gene=goi,
             tx.width=width(tx),
             rna=rna.rpkm,
             pol2.0=pol2.rpkm.0,
             pol2.500=pol2.rpkm.500,
             chroPlus.0=chroPlus.rpkm.0,
             chroPlus.500=chroPlus.rpkm.500,
             chroMinus.0=chroMinus.rpkm.0,
             chroMinus.500=chroPlus.rpkm.500
             )

} # rpkm.4way
#----------------------------------------------------------------------------------------------------
test_rpkm.4way <- function()
{
   goi <- geneSymbols[11]
   tbl <- rpkm.4way(goi)

} # test_rpkm.4way
#----------------------------------------------------------------------------------------------------
bigRun <- function(gois)
{
    tbls <- list()
    for(goi in gois){
       tbl <- tryCatch({
          rpkm.4way(goi)
          },
       error = function(e){
           message(sprintf(" failed with %s", goi))
           print(e)
           data.frame(gene=goi,
                      tx.width=0,
                      rna=0,
                      pol2.0=0,
                      pol2.500=0,
                      chroPlus.0=0,
                      chroPlus.500=0,
                      chroMinus.0=0,
                      chroMinus.500=0)
          } # error
         ) # tryCatch

      tbls[[goi]] <- tbl
      if(tbl$rna > 0)
         print(tbl)
       progress <- grep(goi, gois)[1]
       if((progress %% 100) == 0){
          filename <- sprintf("tbls-%d.RData", progress)
          save(tbls, file=filename)
          }
      } # for

    filename <- sprintf("tbls-all-%s.RData",
                        timestamp <- sub(" ", "", tolower(format(Sys.time(), "%Y.%b.%e-%H:%M"))))
    save(tbls, file=filename)
    browser()
    xyz <- 99

} # bigRun
#----------------------------------------------------------------------------------------------------
if(!interactive()){
    bigRun(geneSymbols)
    }
