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
        #' @param viz logical display tracks on igvR
        #'
        #' @return the rpkm score
        #'
        #' @export
        #'
        rpkm = function(geneSymbol, score.threshold=0.0, viz=FALSE){

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
            exons <- exonsBy(txdb, by="tx")[[tx$tx_id]]
            exons.total.width <- sum(width(exons))

            if(viz){
                tbl.track <- data.frame(chrom=as.character(seqnames(tx)),
                                        start=start(tx),
                                        end=end(tx))
                displayTrack(igv, DataFrameAnnotationTrack("tx", tbl.track, color="random"))
                tbl.track <- data.frame(chrom=as.character(seqnames(exons)),
                                        start=start(exons),
                                        end=end(exons))
                displayTrack(igv, DataFrameAnnotationTrack("exons", tbl.track))
                }

            #------------------------------------------------------------
            # count the reads across this transcript
            #------------------------------------------------------------

            suppressWarnings(
               gr.exons <- import(private$rna.bigwigFile, which=exons, format="bigwig")
               )

            gr.exons.trimmed <- gr.exons[gr.exons$score >= score.threshold]
            total.reads.these.exons <- round(sum(gr.exons.trimmed$score))

            if(viz){
                track <- GRangesQuantitativeTrack("rna-1", gr.exons.trimmed, autoscale=TRUE, min=0.9, max=300, color="random")
                displayTrack(igv, track)
                }

            k <- sum(width(exons))/1000
            rpk <- total.reads.these.exons/k

            if(is.null(private$total.rna.reads)){
               gr.total <- import(private$rna.bigwigFile, format="bigwig")
               private$total.rna.reads <- length(gr.total[gr.total$score >= score.threshold])
               }
            millions.of.total.reads <- private$total.rna.reads/1e6
            rpk / millions.of.total.reads
            } # rpkm

       ) # public

    ) # class
#--------------------------------------------------------------------------------
