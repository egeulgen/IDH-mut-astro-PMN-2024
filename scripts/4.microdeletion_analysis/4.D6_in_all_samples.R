##### Script purpose: Determine deletion status for the 19q13.43 peak identified by GISTIC in all samples
##### Author: Ege Ulgen
##### Date: Dec 2023

library(DESeq2)
library(ggpubr)
library(dplyr)

# identify 19q13.43 deletion status per each sample -----------------------
all_segments <- readRDS("data/selected_data/CN_segments.RDS")

GISTIC_peaks <- read.delim("output/PMN_neg_analysis/all_GISTIC_res/microdel19q/all_lesions.conf_90.txt")
GISTIC_peaks$X <- NULL
GISTIC_peaks <- GISTIC_peaks[!grepl("CN values", GISTIC_peaks$Unique.Name), ]

peak_info <- GISTIC_peaks[, 1:9]

del_19q13.43 <- peak_info[peak_info$Unique.Name == "Deletion Peak 6", ]


scna_gr <- makeGRangesFromDataFrame(all_segments,
                                    seqnames.field = "Chromosome",
                                    start.field = "Start",
                                    end.field = "End",
                                    keep.extra.columns = TRUE)
seqlevelsStyle(scna_gr) <- "UCSC"

del_19q13.43$Wide.Peak.Limits
del_19q13.43$Peak.Limits
microdel_gr <- makeGRangesFromDataFrame(data.frame(chromosome = "chr19",
                                              chromosome_start = 54950535,
                                              chromosome_end = 56859293))

hits <- findOverlaps(scna_gr,
                     microdel_gr,
                     type="any", select="all")

ranges <- microdel_gr[subjectHits(hits)]
mcols(ranges) <- c(mcols(ranges), mcols(scna_gr[queryHits(hits)]))

microdel_df <- data.frame(chr = as.vector(seqnames(ranges)),
                     as.data.frame(mcols(ranges)))
microdel_df <- as.data.frame(microdel_df %>%
                            mutate(segment_start = as.integer(start(ranges(scna_gr[queryHits(hits)])))) %>%
                            mutate(segment_end = as.integer(end(ranges(scna_gr[queryHits(hits)])))) %>%
                            mutate(segment_length_Mb = round((as.numeric((segment_end - segment_start + 1) / 1e6)), digits = 4)) %>%
                            mutate(microdel_start = start(ranges)) %>%
                            mutate(microdel_end = end(ranges)) %>%
                            mutate(chrom = as.character(seqnames(ranges))) %>%
                            rowwise() %>%
                            mutate(microdel_overlap_percent =
                                       round(as.numeric((min(microdel_end, segment_end) - max(segment_start, microdel_start)) / (microdel_end - microdel_start)) * 100, digits = 2)))
# Filter for >= 25% overlap
microdel_df <- microdel_df[microdel_df$microdel_overlap_percent >= 25, ]

all_samples <- unique(microdel_df$patient_barcode)
final_microdel_df <- c()
for (sample_id in all_samples) {
    tmp <- microdel_df[microdel_df$patient_barcode == sample_id, c("patient_barcode", "Segment_Mean")]
    
    if (nrow(tmp) == 1) {
        final_microdel_df <- rbind(final_microdel_df, tmp)
    } else {
        tmp2 <- tmp[1, ]
        tmp2$Segment_Mean <- mean(tmp2$Segment_Mean)
        final_microdel_df <- rbind(final_microdel_df, tmp2)
    }
}

metadata_df <- readRDS("data/selected_data/meta.RDS")
metadata_df$microdel_19q13.43 <- metadata_df$patient %in% final_microdel_df$patient_barcode[final_microdel_df$Segment_Mean < -.3]

table(metadata_df$microdel_19q13.43, metadata_df$class)

saveRDS(metadata_df, "data/selected_data/meta.RDS")

# visulize ----------------------------------------------------------------
tbl <- table(metadata_df$class, metadata_df$microdel_19q13.43)

perc_df <- round(tbl / rowSums(tbl), 4)
perc_df <- as.data.frame(perc_df)

perc_df <- perc_df[perc_df$Var2 == TRUE, ]

p <- ggplot(perc_df, aes(y = Freq, x = Var1, fill = Var1))
p <- p + geom_bar(stat = "identity", position = "stack", aes(color = Var1))
p <- p + geom_text(aes(label = paste0(Freq * 100, "%"), fontface = 2),
                   position = position_stack(vjust = 0.5), size = 4, color = "white")
p <- p + ggsci::scale_fill_lancet()
p <- p + ggsci::scale_color_lancet()
p <- p + ylab("% 19q13.43 microdeletion") + xlab("Class")
p <- p + scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                            expand = c(0, 0))
p <- p + theme_minimal()
p <- p + theme(legend.position = "none",
               axis.text = element_text(size = 11),
               axis.title.y = element_text(face = "bold", size = 14),
               axis.title.x = element_blank())
p

ggsave("output/24.19q13.43_microdel_status.pdf", p, width = 5, height = 7)
