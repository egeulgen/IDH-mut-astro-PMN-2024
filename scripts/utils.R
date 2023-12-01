#' Determine Chromosomal-arm-level SCNA ratio
#' 
#' Calculates and returns the average (weighted arithmetic mean) SCNA ratio used 
#' to determine arm-level CN status
#' 
#'
#' @param SCNA_df input SCNA segments data frame
#' @param cytoband_df input cytobands data frame
#' @param selected_chr selected chromosomes
#'
#' @return a data frame containing: the 'Case_ID', 'arm', 'avg_ratio' columns
determine_chr_arm_level_CN_ratios <- function(SCNA_df, cytoband_df, selected_chr = paste0("chr", c(1, 19))) {

    all_arm_lengths <- data.frame()
    for (chr in selected_chr) {
        tmp <- cytoband_df[cytoband_df$V1 == chr, ]
        tmp$arm <- sub("\\d+.+", "", tmp$V4)
        tmp_p <- tmp[tmp$arm == "p", ]
        tmp_q <- tmp[tmp$arm == "q", ]
        arm_lengths <- data.frame(chrom = chr, arm = paste0(chr, c("p", "q")), start = c(min(tmp_p$V2), min(tmp_q$V2)), end = c(max(tmp_p$V3),
            max(tmp_q$V3)))
        all_arm_lengths <- rbind(all_arm_lengths, arm_lengths)
    }

    all_cases <- unique(SCNA_df$Sample)
    all_arm_ratios_df <- data.frame()
    for (case_ID in all_cases) {
        segments_df <- SCNA_df[SCNA_df$Sample == case_ID, ]
        arm_ratios_df <- all_arm_lengths
        arm_ratios_df$avg_ratio <- 1
        for (chr in selected_chr) {
            for (arm in paste0(chr, c("p", "q"))) {
                tmp_arm <- arm_ratios_df[arm_ratios_df$arm == arm, ]
                c_chr <- segments_df$chr == chr
                c1 <- segments_df$End <= tmp_arm$end & segments_df$End >= tmp_arm$start
                c2 <- segments_df$Start <= tmp_arm$end & segments_df$Start >= tmp_arm$start
                c3 <- segments_df$Start <= tmp_arm$start & segments_df$End >= tmp_arm$end
                arm_segs <- segments_df[c_chr & (c1 | c2 | c3), ]

                arm_segs <- arm_segs[!is.na(arm_segs$Ratio), ]
                arm_segs <- arm_segs[!is.infinite(arm_segs$Ratio), ]
                if (nrow(arm_segs) != 0) {
                  arm_segs$len <- arm_segs$End - arm_segs$Start + 1

                  weighted_mean_ratio <- weighted.mean(arm_segs$Ratio, arm_segs$len/sum(arm_segs$len))
                  arm_ratios_df$avg_ratio[arm_ratios_df$arm == arm] <- weighted_mean_ratio
                }
            }
        }
        arm_ratios_df$Case_ID <- case_ID
        all_arm_ratios_df <- rbind(all_arm_ratios_df, arm_ratios_df)
    }
    chr_arm_df <- all_arm_ratios_df[, c("Case_ID", "arm", "avg_ratio")]
    return(chr_arm_df)
}

#' Produce Genome SCNA Plot for Selected Patient
#'
#' @param donor TCGA barcode of patient
#' @param SCNA_df SCNA segments data frame
#' @param cytoband_df input cytobands data frame
#'
#' @return Genome SCNA Plot
produce_genome_SCNA_plot <- function(donor, SCNA_df, cytoband_df) {
    tm_SCNA_df <- SCNA_df[SCNA_df$patient_barcode == donor, ]

    chr_list <- paste0("chr", 1:22)
    hg38_chr_cumsum <- c(0)
    for (chr in chr_list[-22]) {
        tmp <- cytoband_df[cytoband_df$V1 == chr, ]
        hg38_chr_cumsum <- c(hg38_chr_cumsum, hg38_chr_cumsum[length(hg38_chr_cumsum)] + max(tmp$V3))
    }
    names(hg38_chr_cumsum) <- chr_list

    # add chr end lines
    cumsum_right <- c(0)
    for (chr in chr_list) {
        tmp <- cytoband_df[cytoband_df$V1 == chr, ]
        cumsum_right <- c(cumsum_right, cumsum_right[length(cumsum_right)] + max(tmp$V3))
    }

    plot(NA, xlim = c(0, max(cumsum_right)), ylim = c(-5, 5), xaxs = "i", xaxt = "n", yaxt = "n", xlab = NA, ylab = NA)
    abline(v = cumsum_right, col = "grey")

    # add centromere lines
    cent_pos <- c()
    for (chr in chr_list) {
        tmp <- cytoband_df[cytoband_df$V1 == chr, ]
        cent_pos <- c(cent_pos, hg38_chr_cumsum[chr] + min(tmp$V2[tmp$V5 == "acen"]))
    }

    abline(v = cent_pos, col = "grey", lty = 2)

    chr_sizes <- c()
    for (chr in chr_list) {
        tmp <- cytoband_df[cytoband_df$V1 == chr, ]
        chr_sizes <- c(chr_sizes, max(tmp$V3))
    }

    # x axis
    axis(1, at = hg38_chr_cumsum + chr_sizes/2, labels = chr_list, las = 2)

    # y axis
    axis(2, at = round(seq(-5, 5, 1), 1), las = 2)

    tm_SCNA_df$chromosome2 <- paste0("chr", tm_SCNA_df$Chromosome)

    for (i in seq_len(nrow(SCNA_df))) {
        seg_col <- ifelse(tm_SCNA_df$chromosome2[i] %in% c("chr1", "chr19"), "darkred", "darkblue")
        cur_lwd <- ifelse(tm_SCNA_df$chromosome2[i] %in% c("chr1", "chr19"), 4, 2)

        lines(c(tm_SCNA_df$Start[i] + hg38_chr_cumsum[tm_SCNA_df$chromosome2[i]], tm_SCNA_df$End[i] + hg38_chr_cumsum[tm_SCNA_df$chromosome2[i]]),
            rep(min(5, max(-5, tm_SCNA_df$Segment_Mean[i])), 2), col = seg_col, lwd = cur_lwd)
    }
    title(donor)
}
