# ==============================================================
# 4_create_annotations.R
# Author: Steven Brooks
# Date: 02/19/2025
# --------------------------------------------------------------
# Description:
# This script extracts genomic annotations from the hg38 reference genome
# using Bioconductor packages. It generates annotations for genes, transcripts,
# exons, promoters, coding sequences (CDS), and untranslated regions (UTRs).
#
# Dependencies:
# - Requires `GenomicRanges`, `rtracklayer`, and `TxDb.Hsapiens.UCSC.hg38.knownGene`.
# - The user must specify a reference directory to store annotation files.
#
# Outputs:
# - Saves genomic annotation `.rds` files in the specified directory.
#
# ==============================================================



# ==============================================================
# FUNCTION DEFINITIONS
# ==============================================================

#' Generate genomic annotations and save them as `.rds` files
#'
#' This function extracts annotations for genes, transcripts, exons,
#' promoters, coding sequences (CDS), and untranslated regions (UTRs)
#' from the hg38 reference genome.
#'
#' @param reference_directory Character. Path to the directory where annotations will be stored.
#' @return Saves `.rds` files for each annotation type.
make_annotations <- function(reference_directory) {
  ########### Gene annotations

  # Load the TxDb object
  txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

  # Retrieve genes
  genes <- genes(txdb)
  mcols(genes)$type <- "gene"

  # Retrieve transcripts
  transcripts <- transcripts(txdb)
  mcols(transcripts)$type <- "transcript"

  # Retrieve exons
  exons <- exons(txdb)
  mcols(exons)$type <- "exon"

  # Retrieve promoters (upstream regions of genes)
  promoters <- promoters(txdb, upstream = 2000, downstream = 200)
  mcols(promoters)$type <- "promoter"


  # Retrieve CDS
  cds <- cds(txdb)
  mcols(cds)$type <- "CDS"

  # Retrieve 3' UTRs
  three_utrs <- unlist(threeUTRsByTranscript(txdb))
  mcols(three_utrs)$type <- "UTR"

  # Retrieve 5' UTRs
  five_utrs <- unlist(fiveUTRsByTranscript(txdb))
  mcols(five_utrs)$type <- "UTR"

  # Combine all features into a single GRanges object
  all_features <- c(genes, transcripts, exons, promoters, cds, three_utrs, five_utrs)

  # Sort the combined GRanges object
  all_features <- sort(all_features)


  genome_gr <- keepStandardChromosomes(all_features, pruning.mode = "coarse")
  mcols(genome_gr) <- mcols(genome_gr[, "type", drop = FALSE])

  ########### CpG Islands and TSS


  session <- browserSession()
  genome(session) <- "hg38"
  island_query <- ucscTableQuery(session, table = "cpgIslandExt")
  cpg_islands <- track(island_query)
  cpg_islands <- keepStandardChromosomes(cpg_islands, pruning.mode = "coarse")



  flank_size <- 2000
  # Upstream flank (before the start of the original range)
  upstream_shore <- GRanges(
    seqnames = seqnames(cpg_islands),
    ranges = IRanges(
      start = start(cpg_islands) - flank_size,
      end = start(cpg_islands) - 1
    ),
    strand = strand(cpg_islands)
  )


  # Downstream flank (after the end of the original range)
  downstream_shore <- GRanges(
    seqnames = seqnames(cpg_islands),
    ranges = IRanges(
      start = end(cpg_islands) + 1,
      end = end(cpg_islands) + flank_size
    ),
    strand = strand(cpg_islands)
  )

  # Define shores
  # shores_gr <- flank(cpg_islands, width = 2000, both = TRUE)
  shores_gr <- c(upstream_shore, downstream_shore)

  # Upstream flank (before the start of the original range)
  upstream_shelf <- GRanges(
    seqnames = seqnames(shores_gr),
    ranges = IRanges(
      start = start(shores_gr) - flank_size,
      end = start(shores_gr) - 1
    ),
    strand = strand(shores_gr)
  )

  # Downstream flank (after the end of the original range)
  downstream_shelf <- GRanges(
    seqnames = seqnames(shores_gr),
    ranges = IRanges(
      start = end(shores_gr) + 1,
      end = end(shores_gr) + flank_size
    ),
    strand = strand(shores_gr)
  )

  shelves_gr <- sort(c(upstream_shelf, downstream_shelf))



  # Find open seas
  open_seas_gr <- setdiff(genome_gr, c(cpg_islands, shores_gr, shelves_gr))



  ### Make sure overlaps follow this tiebreaking rule:
  # If a region is annotated by two different areas at the same time:
  #
  #   1) Prioritize CpG island labels.
  # 2) CpG Shores
  # 3) CpG Shelves


  # Find overlaps with CpG islands first and go down

  # template: X = X \ (X ^ Y)
  #
  # X = shores, Y = island
  #
  # X = shelves, Y = shores ^ island

  shores <- setdiff(
    shores_gr,
    intersect(shores_gr, cpg_islands)
  )

  shelves <- setdiff(
    shelves_gr,
    intersect(shelves_gr, c(shores, cpg_islands))
  )

  sea <- setdiff(
    open_seas_gr,
    intersect(open_seas_gr, c(shores, cpg_islands, shelves_gr))
  )


  # Annotate each region type
  mcols(sea)$type <- "Open_Sea"
  mcols(shores)$type <- "Shore"
  mcols(shelves)$type <- "Shelf"
  mcols(cpg_islands)$type <- "CpG_Island"


  cpg_annotation <- c(cpg_islands, shores, shelves, sea)



  # TSS annotations
  # pull the transcripts only
  tss <- genome_gr[genome_gr$type == "transcript"]
  # Shrink down to just the start site.
  ranges(tss) <- IRanges(start = start(tss), end = start(tss) + 1)

  # promoter ~ 1500 upstream, 500 downstream of TSS
  promoter <- reduce(c(
    tss,
    flank(tss, width = 2000, start = TRUE, both = FALSE),
    flank(tss, width = 500, start = FALSE, both = FALSE)
  ))

  mcols(promoter)$type <- "Promoter"

  cpg_annotation <- c(cpg_annotation, promoter)


  saveRDS(genome_gr, file = file.path(reference_directory, "gene_annotations.rds"))
  saveRDS(cpg_annotation, file = file.path(reference_directory, "cpg_annotation.rds"))
}
