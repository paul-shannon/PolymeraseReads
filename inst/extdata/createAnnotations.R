library(annotatr)
aoi <- c("hg38_genes_firstexons")
gr.firstExons <- build_annotations(genome="hg38", annotations=aoi)
length(gr.firstExons) # 272352
f <- "~/github/pol2-ser2-reads/inst/extdata/gr.firstExons.RData"
save(gr.firstExons, file=f)
# tbl.firstExons <- as.data.frame(range(subset(gr.firstExons, symbol=="GATA2")))
# colnames(tbl.firstExons)[1] <- "chrom"
# tbl.firstExons$chrom <- as.character(tbl.firstExons$chrom)
# tbl.firstExons
#   chrom     start       end width strand
# 1  chr3 128481819 128493201 11383      -

