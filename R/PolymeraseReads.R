#' @title PolymeraseReads
#' @description calculates rpkm of pol2-ser2 binding from bigwig data, omitting first exons
#' @name PolymeraseReads
#' @import BiocGenerics
#' @import org.Hs.eg.db
#' @import annotatr
#' @import GenomicFeatures
#' @import TxDb.Hsapiens.UCSC.hg38.knownGene
#' @import rtracklayer
#'
#' @examples
#'   psr <- PolymeraseReads$new(geneSymbols="GATA2")
#'
#' @export

PolymeraseReads = R6Class("PolymeraseReads",

    #--------------------------------------------------------------------------------
    private = list(geneSymbol=NULL,
                   geneID=NULL,
                   txdb=NULL,
                   fileList=c(),
                   total.genome.reads=list(),
                   gr.transcript.hg38=NULL,
                   gr.transcript.hg19=NULL,
                   gr.reads=NULL
                   ),

    #--------------------------------------------------------------------------------
    public = list(
         #' @description
         #' Creates a new instance of this [R6][R6::R6Class] class.
         #' @param fileList named character vector, a list of bigwig files where
         #'   assays (names) can be one or more of rna, pol2, chroPlus, chroMinus
         #' @return a new instance of PolymeraseReads
        initialize = function(fileList){
            assays <- sort(unique(names(fileList)))
            stopifnot(all(assays %in% c("rna", "pol2", "chroPlus", "chroMinus")))
            private$fileList <- fileList
            private$total.genome.reads <- list(rna=NULL, pol2=NULL,
                                               chroPlus=NULL, chroMinus=NULL)
            private$txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
            },
        #------------------------------------------------------------
        #' @description accessor for the object's geneSybmol field
        #' @param geneSymbol character, the HUGO symbol of the gene if interest
        setGeneSymbol = function(geneSymbol){
            private$geneSymbol <- geneSymbol
            private$geneID <- mget(geneSymbol, org.Hs.egSYMBOL2EG, ifnotfound=NA)[[1]]
            },
        #------------------------------------------------------------
        #' @description accessor for the object's geneSybmol field
        #' @return the geneSymbol
        getGeneSymbol = function(){
            private$geneSymbol
            },
        #------------------------------------------------------------
        #' @description accessor for the object's entrez geneID
        #' @return the geneID
        getGeneID = function(){
           private$geneID
           },
        #------------------------------------------------------------
        #' @description find the coordinates of the longest transcript
        #'   in both hg38 and hg19
        #' @param igv igvR instance, default NA
        #' @return a GRanges object
        findLongestTranscript = function(igv=NA){
           if(is.null(private$geneSymbol)) stop()
           if(is.null(private$geneID)) stop()
           tx.by.gene <- transcriptsBy(private$txdb, by="gene")[[private$geneID]]
           if(is.null(tx.by.gene))
               return(NA)
           longest.transcript <- which(width(tx.by.gene) == max(width(tx.by.gene)))[1]
           tx.hg38 <- tx.by.gene[longest.transcript]
           private$gr.transcript.hg19 <- self$lift.hg38.to.hg19(tx.hg38)
           private$gr.transcript.hg38 <- tx.hg38
           if(!is.na(igv)){
               track <- GRangesAnnotationTrack("tx" , tx.hg38, color="random")
               displayTrack(igv, track)
               } # igv
           tx.hg38
           },

        #------------------------------------------------------------
        #' @description return the coordinates of the longest transcript
        #'   in both hg38 and hg19
        #' @return a list of 2 GRanges objects
        getTranscriptCoordinates = function(){
            list(hg38=private$gr.transcript.hg38,
                 hg19=private$gr.transcript.hg19)
            },

        #------------------------------------------------------------
        #' @description liftover GRanges from hg19 to hg38
        #' @param gr  GRanges, the hg19 data structure
        #' @return GRanges
        lift.hg19.to.hg38 = function(gr){
           chain.file <- "hg19ToHg38.over.chain"
           gz.file <- sprintf("%s.gz", chain.file)
           if(!file.exists(chain.file)){
              url <- sprintf("http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/%s", gz.file)
              system(sprintf("curl -O %s", url))
              system(sprintf("gunzip %s", gz.file))
              }
           chain <- import.chain(chain.file)
           x <- liftOver(gr, chain)
           gr.hg38 <- unlist(x)
           seqinfo(gr.hg38) <- SeqinfoForUCSCGenome("hg38")[seqlevels(gr.hg38)]
           gr.hg38
           },

        #' @description liftover GRanges from hg38 to hg18
        #' @param gr  GRanges, the hg38 data structure
        #' @return GRanges
        lift.hg38.to.hg19 = function(gr){
           chain.file <- "hg38ToHg19.over.chain"
           gz.file <- sprintf("%s.gz", chain.file)
           if(!file.exists(chain.file)){
              url <- sprintf("http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/%s", gz.file)
              system(sprintf("curl -O %s", url))
              system(sprintf("gunzip %s", gz.file))
              }
           chain <- import.chain(chain.file)
           x <- liftOver(gr, chain)
           gr.hg19 <- unlist(x)
           seqinfo(gr.hg19) <- SeqinfoForUCSCGenome("hg19")[seqlevels(gr.hg19)]
           gr.hg19
           },

#----------------------------------------------------------------------------------------------------

        #' @description the GRanges of scored and aligned reads of
        #' the longest transcript of geneSymbol
        #' @param assay character, the name of a single assay
        #' @param igv igvR instance, default NA
        #' @param use.igv.roi logical query igv for current region of interest
        #' @return GRanges
        #------------------------------------------------------------
        getReads = function(assay, igv=NA, use.igv.roi=FALSE){
           stopifnot(!is.null(private$gr.transcript.hg38))
           stopifnot(assay %in% names(private$fileList))
           roi <- private$gr.transcript.hg38
           if(!is.na(igv))
               if(use.igv.roi){
                 roi.igv <- getGenomicRegion(igv)
                 roi <- with(roi.igv, GRanges(seqnames=chrom, IRanges(start, end)))
               }
           suppressWarnings({
              gr.reads <- import(private$fileList[assay],
                                 which=roi, format="bigwig")
              })

           if(!is.na(igv)){
               track <- GRangesQuantitativeTrack(assay, gr.reads,
                                                 autoscale=TRUE, color="random")
               displayTrack(igv, track)
               } # igv
           private$gr.reads <- gr.reads
           gr.reads
           },

        #'
        #' @param geneSymbol character HUGO gene symbol name
        #' @param score.threshold numeric ignore reads scoring less than this value
        #' @param assay character either "rna" or "pol2"
        #' @param start.site.avoidance numeric default 500bp, for pol2 only, permits
        #'        rpkm calculation to avoid distortion by
        #' @param igv default NA, a reference to an igvR instance
        #'
        #' @return the rpkm score
        #'
        #' @export
        #'
        rpkm = function(geneSymbol, score.threshold=0.0, assay,
                        start.site.avoidance = 500, igv=NA){

            stopifnot(assay %in% c("rna", "pol2", "chroPlus", "chroMinus"))
            self$setGeneSymbol(geneSymbol)
            tx <- self$findLongestTranscript(igv)
            bigwigFile <- fileList[assay]

            if(assay == "rna"){  # no start.site.avoidance makes sense here
               reads.region <- exonsBy(private$txdb, by="tx")[[tx$tx_id]]
               if(is.null(private$total.genome.reads[[assay]])){
                  message(sprintf("counting total reads in rna file"))
                  gr.total <- import(bigwigFile, format="bigwig")
                  private$total.genome.reads[[assay]] <-
                      length(gr.total[gr.total$score >= score.threshold])
                  }
               } # rna

            if(assay %in% c("pol2", "chroPlus", "chroMinus")){
                  # see chu 2018, chromatin run-on and sequencing, on avoiding
                  # pol2 tss-proximal pileup
               if(as.character(strand(tx)) == "-"){
                  end(tx) <- end(tx) - start.site.avoidance
               } else {
                   start(tx) <- start(tx) + start.site.avoidance
                   }

               reads.region <- tx
               if(is.null(private$total.genome.reads[[assay]])){
                  message(sprintf("counting total reads in %s file", assay))
                  gr.total <- import(bigwigFile, format="bigwig")
                  private$total.genome.reads[[assay]] <-
                          length(gr.total[abs(gr.total$score) >= score.threshold])
                  }
               } # pol2

            suppressWarnings(
               gr.reads <- import(bigwigFile, which=reads.region, format="bigwig")
               )

            gr.reads.trimmed <- gr.reads[abs(gr.reads$score) >= score.threshold]

            if(!is.na(igv)){
               title <- assay
               if(assay != "rna")
                  title <- sprintf("%s.%d", assay, start.site.avoidance)
               track <- GRangesQuantitativeTrack(title, gr.reads.trimmed,
                                                 autoscale=TRUE, color="random")
               displayTrack(igv, track)
               }

            total.reads.in.region <- round(sum(gr.reads.trimmed$score))
            total.reads.in.file <- private$total.genome.reads[[assay]]

            k <- sum(width(reads.region))/1000
            rpk <- total.reads.in.region/k

            millions.of.total.reads <- total.reads.in.file/1e6
            abs(rpk / millions.of.total.reads)
            } # rpkm

       ) # public

    ) # class
#--------------------------------------------------------------------------------
