library(PolymeraseReads)
library(RUnit)

#library(TxDb.Hsapiens.UCSC.hg38.knownGene)
#library(org.Hs.eg.db)

tbl.chr1 <- get(load("~/github/PolymeraseReads/inst/extdata/tbl-rpkm-chr1.RData"))
dim(tbl.chr1) # 3808 10
head(tbl.chr1)
assay.cols <- c("rna", "rna", "pol2.0", "pol2.500", "chroPlus.0", "chroPlus.500",
                "chroMinus.0", "chroMinus.500")
for(coi in assay.cols){
    tbl.chr1[, coi] <- round(tbl.chr1[, coi], digits=2)
    }

head(tbl.chr1)
#      gene strand tx.width   rna pol2.0 pol2.500 chroPlus.0 chroPlus.500 chroMinus.0 chroMinus.500
# 1 A3GALT2      -    14333  0.00   0.89     0.76       3.01         3.11        5.47          5.33
# 2 AADACL3      +    12651  0.00   0.01     0.01       0.00         0.00        0.01          0.01
# 3 AADACL4      +    22992  0.00   0.03     0.03       0.06         0.06        0.05          0.05
# 4   ABCA4      -   128315  0.00   0.13     0.13       0.60         0.60        0.48          0.48
# 5  ABCB10      -    42126 52.90   0.83     0.70       0.38         0.31        2.62          2.42
# 6   ABCD3      +   100275 28.83   0.09     0.06       2.52         1.53        0.03          0.01

tbl.chr1$chroBoth.0 <- with(tbl.chr1, chroPlus.0, chroMinus.0)
tbl.chr1$chroBoth.500 <- with(tbl.chr1, chroPlus.500, chroMinus.500)

geneSymbols <- sort(unique(tbl.chr1$hgnc_symbol))
deleters <- which(nchar(geneSymbols) == 0)
length(deleters) # 0
if(length(deleters) > 0)
   geneSymbols <- geneSymbols[-deleters]
#----------------------------------------------------------------------------------------------------
full.chromosome.pol2.chroBoth.correlation <- function()
{
    cor(tbl.chr1$pol2.0, tbl.chr1$chroBoth.0, method="spearman") # [1] 0.8777124
    with(tbl.chr1, plot(chroBoth.0, pol2.0))
    fivenum(tbl.chr1$pol2.0)  #  0.00  0.00  0.02  0.35 26.99
    fivenum(tbl.chr1$chroBoth.0) # [1]   0.000   0.000   0.020   0.975 927.800

    tbl.tmp <- subset(tbl.chr1, pol2.0 < 15 & chroBoth.0 < 100)
    with(tbl.tmp, cor(pol2.0, chroBoth.0, method="spearman"))  # 0.88
    with(tbl.tmp, cor(pol2.0, chroBoth.0, method="pearson"))  # 0.47
    with(tbl.tmp, plot(chroBoth.0, pol2.0))

    threshold <- 0.0
    tbl.chr1$pol2.present <- tbl.chr1$pol2.0 > threshold
    tbl.chr1$chroBoth.present <- tbl.chr1$chroBoth.0 > threshold
    with(tbl.chr1, cor(chroBoth.present, pol2.present, method="spearman")) # 0.84
    with(tbl.chr1, cor(chroBoth.present, pol2.present, method="pearson"))  # 0.84
    table(tbl.chr1$chroBoth.present)   # FALSE  TRUE
                                       #  1695  2113
    table(tbl.chr1$pol2.present)       # FALSE  TRUE
                                        #  1601  2207

       #----------------------------------------------
       # do this with just chr1:196697434-201020762
       #----------------------------------------------

    #tbl.tmp <- subset(tbl.chr1,

} # full.chromosome.pol2.chroBoth.correlation
#----------------------------------------------------------------------------------------------------
find.cor <- function(tbl.x, colName.1, colName.2, method="pearson")
{
    cor(tbl.x[, colName.1], tbl.x[, colName.2], method=method)

} # find.cor
#----------------------------------------------------------------------------------------------------
test_find.cor <- function()
{
   toi <- tbl.chr1
   dim(toi)  # 3808 12
   round(find.cor(toi, "pol2.0",   "chroBoth.0"), digits=2)    #  0.24
   round(find.cor(toi, "pol2.500", "chroBoth.500"), digits=2)  #  0.32
   round(find.cor(toi, "pol2.500", "chroPlus.500"), digits=2)  #  0.32
   round(find.cor(toi, "pol2.500", "chroMinus.500"), digits=2) #  0.20

   fivenum(subset(tbl.chr1, rna > 0)$rna) #    0.01    0.82    5.08   14.16 1331.86

   toi <- subset(tbl.chr1, rna > 0.5)
   dim(toi)  # 1194 12

   round(find.cor(toi, "pol2.0",   "chroBoth.0"), digits=2)    # 0.28
   round(find.cor(toi, "pol2.500", "chroBoth.500"), digits=2)  # 0.43
   round(find.cor(toi, "pol2.500", "chroPlus.500"), digits=2)  # 0.43
   round(find.cor(toi, "pol2.500", "chroMinus.500"), digits=2) # 0.21

       #----------------------------------
       # rna > 0.5 and tx > 1000
       #----------------------------------
   toi <- subset(tbl.chr1, rna > 0.5 & tx.width > 1000)
   dim(toi)
   round(find.cor(toi, "pol2.0",   "chroBoth.0"), digits=2)    # 0.28
   round(find.cor(toi, "pol2.500", "chroBoth.500"), digits=2)  # 0.44
   round(find.cor(toi, "pol2.500", "chroPlus.500"), digits=2)  # 0.44
   round(find.cor(toi, "pol2.500", "chroMinus.500"), digits=2) # 0.21

       #----------------------------------
       # rna > 0.5 and tx < 1000
       #----------------------------------
   toi <- subset(tbl.chr1, rna > 0.5 & tx.width <= 1000)
   dim(toi)
   round(find.cor(toi, "pol2.0",   "chroBoth.0"), digits=2)    # 0.89
   round(find.cor(toi, "pol2.500", "chroBoth.500"), digits=2)  # 0.24
   round(find.cor(toi, "pol2.500", "chroPlus.500"), digits=2)  # 0.24
   round(find.cor(toi, "pol2.500", "chroMinus.500"), digits=2) # 0.30

       #----------------------------------
       # rna > 5.0 and tx > 1000
       #----------------------------------
   toi <- subset(tbl.chr1, rna > 5.0 & tx.width > 1000)
   dim(toi)  # 757
   round(find.cor(toi, "pol2.0",   "chroBoth.0"), digits=2)    # 0.25
   round(find.cor(toi, "pol2.500", "chroBoth.500"), digits=2)  # 0.39
   round(find.cor(toi, "pol2.500", "chroPlus.500"), digits=2)  # 0.39
   round(find.cor(toi, "pol2.500", "chroMinus.500"), digits=2) # 0.33

       #----------------------------------
       # rna > 5.0 and tx > 10000
       #----------------------------------
   toi <- subset(tbl.chr1, rna > 5.0 & tx.width > 10000)
   dim(toi)
   round(find.cor(toi, "pol2.0",   "chroBoth.0"), digits=2)    # 0.32
   round(find.cor(toi, "pol2.500", "chroBoth.500"), digits=2)  # 0.34
   round(find.cor(toi, "pol2.500", "chroPlus.500"), digits=2)  # 0.34
   round(find.cor(toi, "pol2.500", "chroMinus.500"), digits=2) # 0.25



   as1 <- grep("AS1", tbl.chr1$gene)  # 173
   linc <- grep("LINC", tbl.chr1$gene) # 229
   hyphens <- grep("-", tbl.chr1$gene) # 495

   drops <- unique(c(as1))
   #drops <- unique(c(as1, linc, hyphens))
   length(drops)   # 173

   toi <- tbl.chr1[-drops,]
   toi <- subset(toi, rna > 0.5 & tx.width > 1000)
   dim(toi)   # 1099
   length(as1)

   round(find.cor(toi, "pol2.0",   "chroBoth.0"), digits=2)    # 0.27
   round(find.cor(toi, "pol2.500", "chroBoth.500"), digits=2)  # 0.40
   round(find.cor(toi, "pol2.500", "chroPlus.500"), digits=2)  # 0.40
   round(find.cor(toi, "pol2.500", "chroMinus.500"), digits=2) # 0.22


} # test_find.cor
#----------------------------------------------------------------------------------------------------



