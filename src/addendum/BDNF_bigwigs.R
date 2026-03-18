library(GenomicRanges)
library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg38)

out_dir <- "addendum/data"
if (!dir.exists(out_dir)) dir.create(out_dir)


# EM-seq
em_1 <- readRDS("data/em-seq/em_51.rds")

# Create the GRanges object
# Since it's single-base data, we set end equal to start
gr <- GRanges(
  seqnames = em_1$chr,
  ranges = IRanges(start = em_1$start, end = em_1$start)
)

# Assign the beta value as the "score" track for the BigWig
gr$score <- em_1$beta

# Sort the object (BigWig files require sorted coordinates)
gr <- sort(gr)

# Fetch chromosome lengths for the chromosomes present in your data
seqlengths(gr) <- seqlengths(Hsapiens)[seqlevels(gr)]

export.bw(gr, con = file.path(out_dir, "em_51.bw"))

# Nanopore

nano_1 <- readRDS("data/nanopore/nano_51.rds")


# nano_comb
gr_n <- GRanges(
  seqnames = nano_1$chr,
  ranges = IRanges(start = nano_1$start+1, end = nano_1$start+1)
)

# Assign the beta value as the "score" track for the BigWig
gr_n$score <- (nano_1$NMod + nano_1$NOtherMod) / (nano_1$NMod + nano_1$NOtherMod + nano_1$NCanonical)

# Sort the object (BigWig files require sorted coordinates)
gr_n <- sort(gr_n)

# Fetch chromosome lengths for the chromosomes present in your data
seqlengths(gr_n) <- seqlengths(Hsapiens)[seqlevels(gr_n)]

export.bw(gr_n, con = file.path(out_dir, "nano_51_comb.bw"))


# 5mC

# nano_comb
gr_5mC <- GRanges(
  seqnames = nano_1$chr,
  ranges = IRanges(start = nano_1$start+1, end = nano_1$start+1)
)

# Assign the beta value as the "score" track for the BigWig
gr_5mC$score <- (nano_1$NMod) / (nano_1$NMod + nano_1$NOtherMod + nano_1$NCanonical)

# Sort the object (BigWig files require sorted coordinates)
gr_5mC <- sort(gr_5mC)

# Fetch chromosome lengths for the chromosomes present in your data
seqlengths(gr_5mC) <- seqlengths(Hsapiens)[seqlevels(gr_5mC)]

export.bw(gr_5mC, con = file.path(out_dir, "nano_51_5mC.bw"))



# 5hmC

# nano_comb
gr_5hmC <- GRanges(
  seqnames = nano_1$chr,
  ranges = IRanges(start = nano_1$start+1, end = nano_1$start+1)
)

# Assign the beta value as the "score" track for the BigWig
gr_5hmC$score <- (nano_1$NOtherMod) / (nano_1$NMod + nano_1$NOtherMod + nano_1$NCanonical)

# Sort the object (BigWig files require sorted coordinates)
gr_5hmC <- sort(gr_5hmC)

# Fetch chromosome lengths for the chromosomes present in your data
seqlengths(gr_5hmC) <- seqlengths(Hsapiens)[seqlevels(gr_5hmC)]

export.bw(gr_5hmC, con = file.path(out_dir, "nano_51_5hmC.bw"))