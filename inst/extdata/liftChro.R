library(rtracklayer)
fm <- "GSM3309956_5587_5598_24205_HGC2FBGXX_J_CHR_TGACCA_R1_minus.bw"
fp <- "GSM3309956_5587_5598_24205_HGC2FBGXX_J_CHR_TGACCA_R1_plus.bw"

gr.total <- import(fm, format="bigwig")
chain.file <- "hg19ToHg38.over.chain"
gz.file <- sprintf("%s.gz", chain.file)
if(!file.exists(chain.file)){
   url <- sprintf("http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/%s", gz.file)
   system(sprintf("curl -O %s", url))
   system(sprintf("gunzip %s", gz.file))
   }
chain <- import.chain(chain.file)
x <- liftOver(gr.total, chain)
gr.hg38 <- unlist(x)
seqinfo(gr.hg38) <- SeqinfoForUCSCGenome("hg38")[seqlevels(gr.hg38)]
dups <- findOverlaps(gr.hg38, drop.self=TRUE);
gr.hg38.fixed <- gr.hg38[-queryHits(dups)]
export(gr.hg38.fixed, con="chroMinus-hg38.bw", format="bigwig")



gr.total <- import(fp, format="bigwig")
x <- liftOver(gr.total, chain)
gr.hg38 <- unlist(x)
seqinfo(gr.hg38) <- SeqinfoForUCSCGenome("hg38")[seqlevels(gr.hg38)]
dups <- findOverlaps(gr.hg38, drop.self=TRUE);
gr.hg38.fixed <- gr.hg38[-queryHits(dups)]
export(gr.hg38.fixed, con="chroPlus-hg38.bw", format="bigwig")
