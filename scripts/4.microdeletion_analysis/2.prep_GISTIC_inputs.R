##### Script purpose: Prepare and save seg files for GISTIC analysis
##### Author: Ege Ulgen
##### Date: Dec 2023


metadata_df <- readRDS("data/selected_data/meta.RDS")
all_segs <- readRDS("data/selected_data/CN_segments.RDS")

outdir <- file.path("output", "PMN_neg_analysis", "all_GISTIC_res")

dir.create(outdir, showWarnings = FALSE)

write_seg_file <- function(sel_class) {
    selected_pts <- metadata_df$patient[metadata_df$class == sel_class]
    
    selected_segs <- all_segs[all_segs$patient_barcode %in% selected_pts, ]
    
    # Sample (tab) Chromosome (tab) Start.bp (tab) End.bp (tab) Num.SNPs (tab) Seg.CN
    seg_df <- selected_segs[, c("patient_barcode", "Chromosome", "Start", "End", "Num_Probes", "Segment_Mean")]
    colnames(seg_df) <- c("Sample", "Chromosome", "Start.bp", "End.bp", "Num.SNPs", "Seg.CN")
    write.table(seg_df, file.path(outdir, paste0(sel_class, ".seg")), 
                sep = "\t", quote = FALSE, row.names = FALSE)
    
}

write_seg_file("PMN-hit")
write_seg_file("microdel19q")
write_seg_file("WT")
