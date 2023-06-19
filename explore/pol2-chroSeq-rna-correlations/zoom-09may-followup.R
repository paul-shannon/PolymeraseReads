tbls <- get(load("tbls-chr1-2023.may.08-12:33.RData"))
tbl <- do.call(rbind, tbls)

#----------------------------------------------------------------------------------------------------
correlate <- function(col.1, col.2, sd=FALSE)
{
    col1 <- tbl[, col.1]
    col2 <- tbl[, col.2]
    if(sd){
        col1 <- col1/sd(col1)
        col2 <- col2/sd(col2)
        return(cor(col1, col2, method="spearman"))
        }
    if (!sd){
        return(cor(col1, col2, method="spearman"))
        }

} # correlate
#----------------------------------------------------------------------------------------------------
test_correlate <- function()
{
   correlate("pol2.0", "pol2.500", sd=TRUE)

} # test_correlate
#----------------------------------------------------------------------------------------------------
identify.poised.transcripts <- function()
{
   delta <- with(tbl, pol2.0 - pol2.500)
   fivenum(delta)
   tbl$pol2.poised <- delta
   subset(tbl, pol2.poised > 16)
   pol2.0.sd <- sd(tbl$pol2.0)
   tbl$pol2.0sd <- tbl$pol2.0/sd(tbl$pol2.0)
   tbl$pol2.500sd <- tbl$pol2.500/sd(tbl$pol2.500)
   tbl$pol2.sd.poised <- with(tbl, pol2.0sd - pol2.500sd)
   fivenum(tbl$pol2.sd.poised)

} # identify.poised.transcripts
#----------------------------------------------------------------------------------------------------
