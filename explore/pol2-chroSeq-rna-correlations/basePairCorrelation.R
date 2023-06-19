library(rtracklayer)
data.dir <- "~/github/PolymeraseReads/inst/extdata"
f <- "NS.1828.002.N704---N506.MB_OHRI_Jurkat_POL2_CI_CPM.bw"
gr.p <- import(file.path(data.dir, f))
round(fivenum(gr.p$score), digits=2)  # [1]  0.00  0.09  0.18  0.50 99.14

fivenum(gr.p$score)
tbl.p <- as.data.frame(gr.p) #

colnames(tbl.p)[1] <- "chrom"
tbl.p$chrom <- as.character(tbl.p$chrom)
dim(tbl.p)
extraneous.chroms <- which(nchar(tbl.p$chrom) > 5)
length(extraneous.chroms)
tbl.p <- tbl.p[-extraneous.chroms,]
dim(tbl.p)
tbl.p <-subset(tbl.p, score > 0)
dim(tbl.p)

f <- "chroPlus-hg38.bw"
gr.cp <- import(file.path(data.dir, f))
fivenum(gr.cp$score)  # [1]    1    1    1    2 6018

f <- "chroMinus-hg38.bw"
gr.cm <- import(file.path(data.dir, f))
fivenum(gr.cm$score)
gr.cm$score <- -1 * gr.cm$score
fivenum(gr.cm$score)  #  1    1    1    2 3587

findOverlaps(head(gr.cp, n=500), head(gr.cm, n=500))
  #  [1]        10         251
  #  [2]       173         495

gr.cp[10]$score   # 1
gr.cm[251]$score  # 3

gr.c <- sort(c(gr.cp, gr.cm))
length(gr.c) # 11274481
tbl.c <- as.data.frame(gr.c)

tbl.p.filtered <- subset(tbl.p, score > 0)

tbl.ov <- as.data.frame(findOverlaps(GRanges(tbl.p.filtered), GRanges(tbl.c)))
dim(tbl.ov) # 9098029       2
nrow(tbl.c) # 11130074
nrow(tbl.p.filtered) # 16741128
nrow(tbl.ov)/nrow(tbl.c)  # 0.8174
nrow(tbl.ov)/nrow(tbl.p.filtered)  # 0.54

 #----------------------------------------------------
 # find overlaps in pol2 vs chroBoth in 4M slide 7
 #----------------------------------------------------

chrom.focus <- "chr1"
start.focus <- 196697434
end.focus   <- 201020762

tbl.p.focus <- subset(tbl.p) #, chrom==chrom.focus)
#                             start >= start.focus &
#                             end <= end.focus)
dim(tbl.p.focus)  # 39230
tbl.c.focus <- subset(tbl.c) # , chrom==chrom.focus)
#                             start >= start.focus &
#                             end <= end.focus)
dim(tbl.c.focus)  # 38600
tbl.ov <- as.data.frame(findOverlaps(GRanges(tbl.p.focus), GRanges(tbl.c.focus)))
nrow(tbl.ov)
nrow(tbl.ov)/nrow(tbl.c.focus)  # 0.8476
nrow(tbl.ov)/nrow(tbl.p.focus)  # 0.8340

# gr.cp[10]
# gr.cm[251]
#
# tbl.cp <- as.data.frame(head(gr.cp, n=500))
# tbl.cm <- as.data.frame(head(gr.cm, n=500))
# tbl.c <- rbind(tbl.cp, tbl.cm)
# colnames(tbl.c)[1] <- "chrom"
# tbl.c$chrom <- as.character(tbl.c$chrom)
#
# dim(tbl.c)
#
# tbl.c <- subset(tbl.c, score > 0)
# dim(tbl.c)
# dups <- which(duplicated(tbl.c[, c("chrom", "start")]))
# length(dups)  # 144407
# tbl.c <- tbl.c[-dups,]
# dim(tbl.c)
# dim(tbl.p)
# fivenum(tbl.c$score)
# fivenum(tbl.p$score)
#
# head(tbl.c)
# dups <- which(duplicated(tbl.c[, c("chrom", "start")]))
# length(dups)
#
# for(dup in dups){
#    chrom.dup <- tbl.c$chrom[dup]
#    start.dup <- tbl.c$start[dup]
#    tbl.dup <- subset(tbl.c, chrom==chrom.dup & start==start.dup)
#    score.sum <- sum(tbl.dup$score)
#    }
