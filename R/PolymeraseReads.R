#' @title PolymeraseReads
#' @description calculates rpkm of pol2-ser2 binding from bigwig data, omitting first exons
#' @name PolymeraseReads
#' @import BiocGenerics
#' @import org.Hs.eg.db
#' @import annotatr
#'
#' @examples
#'   psr <- PolymeraseReads$new(geneSymbols="GATA2")
#'
#' @export

PolymeraseReads = R6Class("PolymeraseReads",

    #--------------------------------------------------------------------------------
    private = list(geneSymbol=NULL,
                   geneID=NULL
                   ),

    #--------------------------------------------------------------------------------
    public = list(
         #' @description
         #' Creates a new instance of this [R6][R6::R6Class] class.
         #' @param geneSymbol character, an indentifier for this object
         #' @return a new instance of PolymeraseReads
        initialize = function(geneSymbol){
            private$geneSymbol <- geneSymbol
            private$geneID <- mget("GATA2", org.Hs.egSYMBOL2EG, ifnotfound=NA)[[1]]
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
        #' @description the firstexons and their range
        #' @return the geneID
        getFirstExons = function(){
           gr.firstExons <- get(load(system.file(package="PolymeraseReads",
                                                 "extdata", "gr.firstExons.RData")))

           gr.sub <- subset(gr.firstExons, symbol==private$geneSymbol)
           tbl.range <- as.data.frame(range(gr.sub))
           colnames(tbl.range)[1] <- "chrom"
           tbl.range$chrom <- as.character(tbl.range$chrom)
           browser()
           firstExonsEnd <- tbl.range$end
           if(tbl.range$strand == "-")
               firstExonsEnd <- tbl.range$start
           tbl.exons <- as.data.frame(gr.sub)
           colnames(tbl.exons)[1] <- "chrom"
           tbl.exons$chrom <- as.character(tbl.exons$chrom)
           list(range=tbl.range, exons=tbl.exons, firstExonsEnd=firstExonsEnd)
           }

       ) # public

    ) # class
#--------------------------------------------------------------------------------
