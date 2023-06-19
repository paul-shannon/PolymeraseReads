library(RUnit)
library(org.Hs.eg.db)
#-------------------------------------------------------------------------------------
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


#-------------------------------------------------------------------------------------
add.zscores <- function(column.title)
{
   stopifnot(column.title %in% colnames(tbl))
   vec <- tbl[, column.title]
   vec.sd <- sd(vec)
   vec.z <- vec/vec.sd
   new.column.title <- sprintf("%s.z", column.title)
   tbl[, new.column.title] <<- vec.z

} # add.zscores
#-------------------------------------------------------------------------------------
test_and.add.zscores <- function()
{
   message(sprintf("--- test_add.zscores"))
   add.zscores("rna")
   checkTrue("rna.z" %in% colnames(tbl))
   checkEqualsNumeric(subset(tbl, rna > 1331.85)$rna.z, 34.62032, tol=1e-3)

   add.zscores("pol2.0")
   add.zscores("pol2.500")
   add.zscores("chroPlus.0")
   add.zscores("chroPlus.500")
   add.zscores("chroMinus.0")
   add.zscores("chroMinus.500")

} # test_add.zscores
#-------------------------------------------------------------------------------------
# median half-life of human mRNA is 10 hours.
# TF mRNA has half-lives of < 2 hours
# pol2 reads on the transcript body are also short lived: it is a dynamic process
# do we see that TFs correlate better with pol2 and chro-seq than everything else?
explore.tf.correlation <- function()
{
  #cor.method <- "spearman"
  cor.method <- "pearson"
  with(subset(tbl, gene %in% tfs.all), cor(pol2.0, rna, method=cor.method))    # 0.270
  with(subset(tbl, ! gene %in% tfs.all), cor(pol2.0, rna, method=cor.method))  # 0.132

  with(subset(tbl, strand=="+" & gene %in% tfs.all),
       cor(chroPlus.500, rna, method=cor.method))    # 0.282

  with(subset(tbl, strand=="+" & ! gene %in% tfs.all),
       cor(chroPlus.500, rna, method=cor.method))    # 0.257

  with(subset(tbl, strand=="-" & gene %in% tfs.all),
       cor(chroMinus.500, rna, method=cor.method))    # -0.006

  with(subset(tbl, strand=="-" & ! gene %in% tfs.all),
       cor(chroMinus.500, rna, method=cor.method))    # 0.121

    # now look more closely
  tbl.filtered <- subset(tbl, gene %in% tfs.all & rna <100 & pol2.500 <6)
  tbl.filtered <- subset(tbl, gene %in% tfs.all & rna < 100)
  with(tbl.filtered, plot(pol2.500, rna, main="TFs only"))
  fit <- lm(rna ~ 0 + pol2.500, data=tbl.filtered)
  abline(fit, col="red")
  summary(fit)

  tbl.filtered <- subset(tbl, !gene %in% tfs.all & rna <100)
  tbl.filtered <- subset(tbl, !gene %in% tfs.all) #  & rna <100 & pol2.500 <6)
  with(tbl.filtered, plot(pol2.500, rna, main="no TFs"))
  fit <- lm(rna ~ 1 + pol2.500, data=tbl.filtered)
  abline(fit, col="red")
  summary(fit)


} # explore.tf.correlation
#-------------------------------------------------------------------------------------
find.mismatched.pol2.rna <- function()
{
  with(tbl, cor(pol2.500.z, rna.z, method="spearman"))  # 0.737

} # find.mismatched.pol2.rna
#-------------------------------------------------------------------------------------

hist(tbl$rna, breaks=seq(0, 1500, by=20))
hist(tbl$pol2.500, breaks=seq(0, 1500, by=20))
hist(tbl$chroPlus.500, breaks=seq(0, 1500, by=20))

fivenum(tbl$rna)           #     0.000000    0.000000    0.000000     2.467815   1331.857143
fivenum(tbl$chroPlus.500)  #     0.000000    0.0000000   0.01325992   0.69179247  372.58606638

round(fivenum(tbl$rna), digits=4)           #   0.0000    0.0000    0.0000    2.4678 1331.8571
round(fivenum(tbl$chroPlus.500), digits=4)  #    0.0000   0.0000    0.0133     0.6918  372.5861

pol2.bot <- fivenum(tbl$pol2.500)[1]
pol2.top <-    fivenum(tbl$pol2.500)[5]
chroPlus.bot <-   fivenum(tbl$chroPlus.500)[1]
chroPlus.top <-   fivenum(tbl$chroPlus.500)[5]
chroMinus.bot <-   fivenum(tbl$chroMinus.500)[1]
chroMinus.top <-   fivenum(tbl$chroPlus.500)[5]

tbl.sub <- subset(tbl, strand=="-" &
                       rna >= 0 &
                       tx.width > 1000 &
                       pol2.500 >= pol2.bot & pol2.500 <= pol2.top &
                       chroMinus.500 >= chroMinus.bot & chroMinus.500 <= chroMinus.top)
dim(tbl.sub)
round(with(tbl.sub, cor(chroMinus.500, pol2.500, method="spearman")), digits=3) # 0.763
round(with(tbl.sub, cor(chroMinus.500, pol2.500, method="pearson")), digits=3) # 0.34


# fivenums      1,4       1,5    3,5     1,5     3,5
# rna           1          1      1       0       0
# size          45        550    520    1298     852
# spearman      0.699     0.57   0.517   0.767   0.585
# pearson       0.459     0.463  0.461   0.316   0.29

round(with(tbl.sub, cor(chroPlus.500, rna, method="spearman")), digits=3) # 0.751
with(subset(tbl.sub, strand=="+"), plot(chroPlus.500, rna))


 > with(subset(tbl, strand=="+" & rna > 5), cor(chroPlus.500, rna, method="spearman"))
 [1] 0.3639371
 > with(subset(tbl, strand=="+" & rna > 1), cor(chroPlus.500, rna, method="spearman"))
 [1] 0.4130893
 > with(subset(tbl, strand=="+"), cor(chroPlus.500, rna, method="spearman"))
 [1] 0.7507425
 > with(subset(tbl, strand=="-"), cor(chroMinus.500, rna, method="spearman"))
 [1] 0.7306893
 > with(subset(tbl, strand=="-"), cor(chroMinus.0, rna, method="spearman"))
 [1] 0.7385342
 > with(subset(tbl, strand=="+"), cor(chroPlus.0, rna, method="spearman"))
 [1] 0.7538335
 > with(subset(tbl, strand=="+"), cor(chroPlus.0, rna, method="pearson"))
 [1] 0.3554539
 > with(subset(tbl, strand=="+"), cor(pol2.500, rna, method="spearman"))
 [1] 0.6141025
 > with(subset(tbl, strand=="-"), cor(pol2.500, rna, method="spearman"))
 [1] 0.6005395
 > with(subset(tbl, strand=="-"), cor(pol2.0, rna, method="spearman"))


