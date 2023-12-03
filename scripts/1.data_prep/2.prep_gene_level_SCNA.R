##### Script purpose: Prepare gene-level SCNA data 
##### Author: Ege Ulgen 
##### Date: Dec 2023

library(Homo.sapiens)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(dplyr)

SCNA_df <- readRDS("data/all_primary_complete/CN_segments.RDS")
SCNA_df <- SCNA_df[, c("patient_barcode", "Chromosome", "Start", "End", "Segment_Mean")]

# gene-level SCNA ---------------------------------------------------------
scna_gr <- makeGRangesFromDataFrame(SCNA_df, seqnames.field = "Chromosome", start.field = "Start", end.field = "End", keep.extra.columns = TRUE)

seqlevelsStyle(scna_gr) <- "UCSC"

genes_gr <- suppressWarnings(genes(TxDb.Hsapiens.UCSC.hg38.knownGene))

### Overlaps
hits <- findOverlaps(scna_gr, genes_gr, type = "any", select = "all")

### Genes df
ranges <- genes_gr[subjectHits(hits)]
mcols(ranges) <- c(mcols(ranges), mcols(scna_gr[queryHits(hits)]))

genes_df <- data.frame(chr = as.vector(seqnames(ranges)), as.data.frame(mcols(ranges)))
genes_df <- genes_df %>%
    mutate(segment_start = as.integer(start(ranges(scna_gr[queryHits(hits)])))) %>%
    mutate(segment_end = as.integer(end(ranges(scna_gr[queryHits(hits)])))) %>%
    mutate(segment_length_Mb = round((as.numeric((segment_end - segment_start + 1)/1e+06)), digits = 4)) %>%
    mutate(transcript_start = start(ranges)) %>%
    mutate(transcript_end = end(ranges)) %>%
    mutate(chrom = as.character(seqnames(ranges))) %>%
    rowwise() %>%
    mutate(
        transcript_overlap_percent = round(
            as.numeric((min(transcript_end, segment_end) - max(segment_start, transcript_start))/(transcript_end - transcript_start)) * 100, 
            digits = 2)
        )
genes_df <- as.data.frame(genes_df)
### Filter for >= 25% overlap
genes_df <- genes_df[genes_df$transcript_overlap_percent >= 25, ]

### add gene symbols
all_syms_tbl <- as.data.frame(org.Hs.egSYMBOL)
genes_df$symbol <- all_syms_tbl$symbol[match(genes_df$gene_id, all_syms_tbl$gene_id)]
genes_df <- genes_df[!is.na(genes_df$symbol), ]

saveRDS(genes_df, "data/all_primary_complete/CN_gene_level.RDS")
