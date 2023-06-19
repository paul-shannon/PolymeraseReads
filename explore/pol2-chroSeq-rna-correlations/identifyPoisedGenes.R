library(RUnit)
library(org.Hs.eg.db)
#-------------------------------------------------------------------------------------
if(!exists("tbl")){
   tbls <- get(load("tbls-chr1-2023.may.08-12:33.RData"))
   tbl <- do.call(rbind, tbls)
   dim(tbl)
   table(tbl$strand)   #         -    +
                       # 1246 1264 1298
  tbl.tfs <- get(load("~/github/MotifDb/inst/extdata/tfs/tfs-1683-lambert.RData"))
  dim(tbl.tfs)   # 1683 5
  length(intersect(tbl.tfs$Gene, tbl$gene)) # 128
  tfs.all <- sort(unique(select(org.Hs.eg.db, keys="GO:0003700", keytype="GOALL", columns="SYMBOL")$SYMBOL))
  length(tfs.all) # 1429
   length(intersect(tfs.all, tbl$gene)) # 109
   }
#-------------------------------------------------------------------------------------
# identify poised transcripts by comparing first 500 bp against remainder
identifyPoisedGenes <- function()
{
   tbl.poised <- subset(tbl, pol2.0 > 1)
   dim(tbl.poised)   # 497 10
     # if the delta is small, most of the reads are in the 5' 500 bp
   round(fivenum(tbl.poised$pol2.0), digits=2)
   hist(tbl.poised$pol2.0)
   tbl.poised$pol2.delta <- with(tbl.poised, pol2.0 - pol2.500)
   round(fivenum(tbl.poised$pol2.delta), digits=2) #  -3.87  0.00  0.00  0.02 16.19

      # poised genes are, in principle, those with nearly identical
      # pol2.0 and pol2.500
   subset(tbl.poised, abs(tbl.poised$pol2.delta) < 0.01)

   with(tbl.poised, hist(pol2.0)) #, pol2.500))
   quartz()
   with(tbl.poised, hist(pol2.500)) #, pol2.500))
   boxplot(tbl.poised$pol2.delta)
   boxplot(tbl.poised$pol2.delta)
   hist(tbl.poised$pol2.delta)

   dim(subset(tbl.poised, pol2.delta < 0))  # 837
   subset(tbl.poised, pol2.delta < -3.8)
   subset(tbl.poised, pol2.delta < -3.8)


} # identifyPoisedGenes
#-------------------------------------------------------------------------------------
