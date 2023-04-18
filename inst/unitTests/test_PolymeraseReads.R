library(RUnit)
library(PolymeraseReads)
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
    test_ctor()

} # runTests
#----------------------------------------------------------------------------------------------------
test_ctor <- function()
{
    message(sprintf("--- test_ctor"))

    geneSymbol <- "GATA2"
    psr <- PolymeraseReads$new(geneSymbol)
    checkTrue(all(c("R6", "PolymeraseReads") %in% class(psr)))
    checkEquals(psr$getGeneSymbol(), geneSymbol)
    checkEquals(psr$getGeneID(), "2624")

} # test_ctor
#----------------------------------------------------------------------------------------------------
test_firstExons <- function()
{
    message(sprintf("--- test_firstExons"))
    geneSymbol <- "GATA2"
    psr <- PolymeraseReads$new(geneSymbol)

    x <- psr$getFirstExons()
    checkTrue(all(c("range", "exons", "firstExonsEnd") %in% names(x)))
    x$range
    x$firstExonsEnd
    if(FALSE){
      igv <- start.igv(geneSymbol, "hg38")
      track <- DataFrameAnnotationTrack("firstExons", x$exons, color="brown")
      displayTrack(igv, track)
      tbl.track <- x$range
      tbl.track$start <- x$firstExonsEnd
      tbl.track$end <- x$firstExonsEnd
      track <- DataFrameAnnotationTrack("firstExonEnd", tbl.track, color="black")
      displayTrack(igv, track)
      }


} # test_firstExons
#----------------------------------------------------------------------------------------------------
if(!interactive())
    runTests()
