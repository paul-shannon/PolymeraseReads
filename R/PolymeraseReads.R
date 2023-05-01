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
                   pol2.bigwigFile=NULL,
                   rna.bigwigFile=NULL,
                   total.rna.reads=NULL,
                   total.pol2.reads=NULL,
                   gr.transcript=NULL,
                   gr.reads=NULL
                   ),

    #--------------------------------------------------------------------------------
    public = list(
         #' @description
         #' Creates a new instance of this [R6][R6::R6Class] class.
         #' @param geneSymbol character, an indentifier for this object
         #' @param pol2.bigwigFile character, full path to the aligned reads file
         #' @param rna.bigwigFile character, full path to the aligned reads file
         #' @return a new instance of PolymeraseReads
        initialize = function(geneSymbol, pol2.bigwigFile, rna.bigwigFile=NA){
            private$pol2.bigwigFile <- pol2.bigwigFile
            private$rna.bigwigFile <- rna.bigwigFile
            private$geneSymbol <- geneSymbol
            private$geneID <- mget(geneSymbol, org.Hs.egSYMBOL2EG, ifnotfound=NA)[[1]]
            if(is.na(private$geneID))
                stop("error.  found no entrez geneID for %s", geneSymbol)
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
        #' @description the coordinates of the longest transcript
        #' @return a list
        getTranscriptCoordinates = function(){
           txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
           tx.by.gene <- transcriptsBy(txdb, by="gene")[[private$geneID]]
           longest.transcript <- which(width(tx.by.gene) == max(width(tx.by.gene)))
           private$gr.transcript <- tx.by.gene[longest.transcript]
           private$gr.transcript
           },

        #' @description the GRanges of scored and aligned reads of
        #' the longest transcript of geneSymbol
        #' @return GRanges
        #------------------------------------------------------------
        getReads = function(){
           stopifnot(!is.null(private$gr.transcript))
           suppressWarnings({
              private$gr.reads <- import(private$pol2.bigwigFile,
                                         which=private$gr.transcript, format="bigwig")
              })

           private$gr.reads
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

            stopifnot(assay %in% c("rna", "pol2"))

            # get the coordinates and width of the longest transcript
            # we presume the longest will be least likely to overcount
            # reads per transcript

            id <- mget(geneSymbol, org.Hs.egSYMBOL2EG, ifnotfound=NA)[[1]]
            stopifnot(!is.na(id))
            txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
            tx.by.gene <- transcriptsBy(txdb, by="gene")[[id]]
            if(is.null(tx.by.gene)){
                message(sprintf("no transcripts for %s", geneSymbol))
                return(0)
                }
            longest.transcript <- which(width(tx.by.gene) == max(width(tx.by.gene)))  # 10
            tx <- tx.by.gene[longest.transcript]

            if(!is.na(igv)){
               tbl.track <- data.frame(chrom=as.character(seqnames(tx)),
                                       start=start(tx),
                                       end=end(tx))
               displayTrack(igv, DataFrameAnnotationTrack("tx", tbl.track, color="random"))
               }
            bigwigFile <- switch(assay,
                                 "rna" = private$rna.bigwigFile,
                                 "pol2" = private$pol2.bigwigFile)
            if(assay == "rna"){
               reads.region <- exonsBy(txdb, by="tx")[[tx$tx_id]]
               if(is.null(private$total.rna.reads)){
                  message(sprintf("counting total reads in rna file"))
                  gr.total <- import(bigwigFile, format="bigwig")
                  private$total.rna.reads <- length(gr.total[gr.total$score >= score.threshold])
                  }
               total.reads.in.file <- private$total.rna.reads
               } # rna

            if(assay == "pol2"){
                  # see chu 2018, chromatin run-on and sequencing, on avoiding
                  # pol2
               if(as.character(strand(tx)) == "-"){
                  end(tx) <- end(tx) - start.site.avoidance
               } else {
                   start(tx) <- start(tx) + start.site.avoidance
                   }

               reads.region <- tx
               if(is.null(private$total.pol2.reads)){
                  message(sprintf("counting total reads in pol2 file"))
                  gr.total <- import(bigwigFile, format="bigwig")
                  private$total.pol2.reads <-
                          length(gr.total[gr.total$score >= score.threshold])
                  }
               total.reads.in.file <- private$total.pol2.reads
               } # pol2

            suppressWarnings(
               gr.reads <- import(bigwigFile, which=reads.region, format="bigwig")
               )

            gr.reads.trimmed <- gr.reads[gr.reads$score >= score.threshold]
            if(!is.na(igv)){
               title <- switch(assay,
                               "pol2" = sprintf("%s %d", assay, start.site.avoidance),
                               "rna"  = "rna")
               track <- GRangesQuantitativeTrack(title, gr.reads.trimmed,
                                                 autoscale=TRUE, color="random")
               displayTrack(igv, track)
               }

            total.reads.in.region <- round(sum(gr.reads.trimmed$score))

            k <- sum(width(reads.region))/1000
            rpk <- total.reads.in.region/k

            millions.of.total.reads <- total.reads.in.file/1e6
            rpk / millions.of.total.reads
            } # rpkm

       ) # public

    ) # class
#--------------------------------------------------------------------------------
