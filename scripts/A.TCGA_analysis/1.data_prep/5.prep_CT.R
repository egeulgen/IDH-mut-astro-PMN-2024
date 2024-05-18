##### Script purpose: Detect chromothripsis events via CTLPScanner
##### Author: Ege Ulgen 
##### Date: Dec 2023

####################################################################################################
####################################################################################################
#
#	Scan statistic to detect clustering of copy number status switch. Three parameters:
#		Minimum segment size
#       Window size
#		Log2 signal distance between adjacent segments as status change
#
#	Haoyang
#	2012-08-01
#
# README ----------------------------------------------------------------
##############################################################################
#	CTLPScanner
#	Chromothripsis-like patterns detection algorithm
#	input: segmented array data
#	output: chromothripsis-like regions
#	© 2015 Cai Laboratory
##############################################################################
#
#
#Please follow the below steps to run the script:
#
#1. This script runs on segmented files. Put the name of files that you want to process in “inputlist???.
#Here are two examples: “GSM681792,segments.tab??? and “GSM681810,segments.tab???.
#
#2. Copy segmented files into the same folder. The file format is as follows:
#
#sample: sample ID
#chro: Chromosome
#start: segment start position
#stop: segment stop position
#mean: log2 value
#probes: number of probes in this segment(optional)
#
#3. Run “CTLP_detection.R???. There are two parameters:
#
#Minimum segment size (default 10Kb): The segment will be ignored if its size is smaller than this value. In our case, we choose 10Kb.
#Log2 value between adjacent segments (default 0.3): This parameter is used to filter out “hypersegmented??? segments. They are more likely to be segmentation noise. If the log2 signal between two adjacent segments is smaller than this value, no status change is considered. We choose 0.3 for Affymetrix SNP arrays.
#
#4. The script will generate a result file (Results.txt).
#In this file, each row represents one Chromosome, including information of the window with the highest likelihood ratio.
#It contains the following columns:
#
#arrayID
#WinNo
#WinSize: window size
#Chrom: Chromosome
#Start: window start position
#End: window end position
#ExpSwitchNo
#SwitchNo: calculated status switch number
#LR : likelihood ratio
#Pvalue
#
#5. Choose the appropriate thresholds to “call chromothripsis-like patterns.
#According to the training set, the default thresholds are: SwitchNo >= 20 and log10(LR) >= 8.
#With the included two examples, Chromosome 2 of GSM681792 and Chromosome 10 of GSM681810 are identified as chromothripsis.
####################################################################################################
####################################################################################################

# prep data ---------------------------------------------------------------
metadata_df <- readRDS("data/selected_data/meta.RDS")

SCNA_df <- readRDS("data/selected_data/CN_segments.RDS")

SCNA_df <- SCNA_df[SCNA_df$Chromosome %in% 1:22, ]
SCNA_df$Chromosome <- as.numeric(SCNA_df$Chromosome)

SCNA_df <- SCNA_df[, c("patient_barcode", "Chromosome", "Start", "End", "Segment_Mean")]

# prep params -------------------------------------------------------------

# Set segment min size (default 10KB)
minSegSize <- 1e4
# Set log2 signal distance between adjacent segments as status change (default 0.3)
segDistance <- 0.3

cytoband_df <- read.delim("data/hg38_cytoBand.txt", header = FALSE)
chromEnd <- c()
for (chr in paste0("chr", 1:22)) {
    tmp <- cytoband_df[cytoband_df$V1 == chr, ]
    chromEnd <- c(chromEnd, max(tmp$V3) + 1)
}

scanMatrix <- matrix(c("patient_barcode", "WinNo", "WinSize", "Chrom", "Start", "End", "ExpSwitchNo", "SwitchNo", "LR", "Pvalue"), nrow = 1)
for (i in 1:22){
    scanMatrix <- rbind(scanMatrix, 
                        c("", nrow(scanMatrix), chromEnd[i] / 1e6, i, 1, chromEnd[i], "","","",""))
}

####################################################################################################
#	Create sliding window
####################################################################################################

windowSize <- sort(c(seq(1e6, 2.5e8, by = 1e7), chromEnd))
windowStep <- 1e6

for (i in 1:length(windowSize)) {
    for (j in 1:22) {
        chromWindowEnd <- 0
        if (windowSize[i] < chromEnd[j]) {
            scanMatrix <- rbind(scanMatrix, c("", nrow(scanMatrix), windowSize[i] / 1e6, j , 1, windowSize[i], "", "", "", ""))
            chromWindowEnd <- windowSize[i] + windowStep
            while (chromWindowEnd < chromEnd[j]){
                scanMatrix <- rbind(scanMatrix, c("", nrow(scanMatrix), windowSize[i] / 1e6, j, chromWindowEnd - windowSize[i] + 1, chromWindowEnd, "", "", "", ""))
                chromWindowEnd <- chromWindowEnd + windowStep
            }
            if ((chromEnd[j] - windowSize[i]) > 0){
                scanMatrix <- rbind(scanMatrix, c("", nrow(scanMatrix), windowSize[i] / 1e6, j, chromEnd[j] - windowSize[i], chromEnd[j], "", "", "", ""))
            }
        }
    }
}
scanMatrix <- rbind(scanMatrix,c("", "ALL", 0, "ALL", "", "", "", "", "", ""))

####################################################################################################
#	Create output folder
####################################################################################################
ctlp_df <- data.frame()

all_patient_barcodes <- unique(SCNA_df$patient_barcode)

for (i in 1:length(all_patient_barcodes)){
    cat("Working on", i, all_patient_barcodes[i], "      \r")
    tmpData <- c()
    arrayMatrix <- scanMatrix
    arraySeg <- SCNA_df[SCNA_df$patient_barcode == all_patient_barcodes[i], ]
    
    tmpData <- arraySeg[which((arraySeg[,2] != "NA") & (arraySeg[,3] != "NA") & (arraySeg[,4] != "NA") & (arraySeg[,5] != "NA")),]
    
    arraySwitchNo <- 0
    for (k in 1:22){
        curSig <- c()
        chromSeg <- tmpData[which((tmpData[,2] == k) & ((tmpData[,4] - tmpData[,3]) >= minSegSize)),]
        if (!is.null(chromSeg)){
            if (nrow(chromSeg) > 2){
                chromSeg <- chromSeg[order(chromSeg[,3]),]
                
                for (segNo in 1:nrow(chromSeg)){
                    if (is.null(curSig)){
                        curSig <- chromSeg[segNo, 5]
                    }
                    else {
                        if (abs(curSig - chromSeg[segNo, 5]) >= segDistance){
                            arraySwitchNo <- arraySwitchNo + 1
                            curSig <- chromSeg[segNo, 5]
                        }
                    }
                }
            }
        }
    }
    arrayMatrix[nrow(arrayMatrix),8] <- arraySwitchNo
    arrayMatrix[nrow(arrayMatrix),7] <- arraySwitchNo
    
    
    for (j in 2:(nrow(arrayMatrix)-1)){
        arrayMatrix[j,1] <- all_patient_barcodes[i]
        chromSeg <- tmpData[which(tmpData[,2] == as.numeric(arrayMatrix[j,4])),]
        
        switchNo <- 0
        curSig <- c()
        
        if (!is.null(chromSeg)){
            if (nrow(chromSeg) > 2){
                chromSeg <- chromSeg[order(chromSeg[,3]),]
                
                winSeg <- chromSeg[which(((chromSeg[,3] >= as.numeric(arrayMatrix[j,5])) & (chromSeg[,3] <= as.numeric(arrayMatrix[j,6]))) | ((chromSeg[,4] >= as.numeric(arrayMatrix[j,5])) & (chromSeg[,4] <= as.numeric(arrayMatrix[j,6])))),]
                winSeg <- winSeg[which((winSeg[,4] - winSeg[,3]) >= minSegSize),]
                
                if (!is.null(winSeg)){
                    if (nrow(winSeg) > 0){
                        for (segNo in 1:nrow(winSeg)){
                            if (is.null(curSig)){
                                curSig <- winSeg[segNo,5]
                            }
                            else {
                                if (abs(curSig - winSeg[segNo,5]) >= segDistance){
                                    switchNo = switchNo + 1
                                    curSig <- winSeg[segNo,5]
                                }
                            }
                        }
                    }
                }
                arrayMatrix[j,7] <- as.numeric(arrayMatrix[j,3])/2868*arraySwitchNo
                arrayMatrix[j,8] <- switchNo
                
                ####################################################################################################
                #	Calculate P value
                ####################################################################################################
                
                simuHigher <- 0
                for (n in 1:1){
                    simuData <- tmpData
                    simuData[,5] <- sample(simuData[,5])
                    simuChromSeg <- simuData[which(simuData[,2] == as.numeric(arrayMatrix[j,4])),]
                    simuChromSeg <- simuChromSeg[order(simuChromSeg[,3]),]
                    simuSwitchNo <- 0
                    simuCurSig <- c()
                    
                    simuWinSeg <- simuChromSeg[which(((simuChromSeg[,3] >= as.numeric(arrayMatrix[j,5])) & (simuChromSeg[,3] <= as.numeric(arrayMatrix[j,6]))) | ((simuChromSeg[,4] >= as.numeric(arrayMatrix[j,5])) & (simuChromSeg[,4] <= as.numeric(arrayMatrix[j,6])))),]
                    simuWinSeg <- simuWinSeg[which((simuWinSeg[,4] - simuWinSeg[,3]) >= minSegSize),]
                    
                    if (!is.null(simuWinSeg)){
                        if (nrow(simuWinSeg) > 0){
                            for (simuSegNo in 1:nrow(simuWinSeg)){
                                if (is.null(simuCurSig)){
                                    simuCurSig <- simuWinSeg[simuSegNo,5]
                                }
                                else {
                                    if (isTRUE(abs(simuCurSig - simuWinSeg[simuSegNo,5]) >= segDistance)){
                                        simuSwitchNo = simuSwitchNo + 1
                                        simuCurSig <- simuWinSeg[simuSegNo,5]
                                    }
                                }
                            }
                        }
                    }
                    if (simuSwitchNo > switchNo){
                        simuHigher <- simuHigher + 1
                    }
                }
                arrayMatrix[j,10] <- simuHigher/1
            }
        }
    }
    
    
    ####################################################################################################
    #	Scan Statistic
    ####################################################################################################
    
    for (j in 2:(nrow(arrayMatrix)-1)){
        LR <- 0
        nz <- as.integer(as.numeric(arrayMatrix[j,8]))
        uz <- as.integer(as.numeric(arrayMatrix[j,7]))
        ng <- as.integer(as.numeric(arrayMatrix[nrow(arrayMatrix),8]))
        ug <- as.integer(as.numeric(arrayMatrix[nrow(arrayMatrix),7]))
        
        if ((!is.na(uz)) & (!is.na(uz)) & (!is.na(ng)) & (!is.na(ug))){
            if (uz > 0){
                if ((nz/uz) > ((ng-nz)/(ug-uz))){
                    LR <- (nz/uz)^nz*((ng-nz)/(ug-uz))^(ng-nz)/(ng/ug)^ng
                }
                else {
                    LR <- 1
                }
            }
        }
        arrayMatrix[j,9] <- LR
    }
    
    
    ####################################################################################################
    #	Find the window with the highest LR. Pick the highest CNA number if the same LR.
    ####################################################################################################
    sortMatrix <- c()
    for (l in 1:22){
        chromWindow <- arrayMatrix[which(as.numeric(arrayMatrix[,4]) == l),]
        sortMatrix <- chromWindow[order(as.numeric(chromWindow[,9]),-as.numeric(chromWindow[,3]),decreasing=T),]
        ctlp_df <- rbind(ctlp_df, 
                         sortMatrix[1,])
    }
}

colnames(ctlp_df) <- c("patient_barcode", "WinNo", "WinSize", "Chrom", "Start", "End", "ExpSwitchNo", "SwitchNo", "LR", "Pvalue")
ctlp_df$WinNo <- as.numeric(ctlp_df$WinNo)
ctlp_df$WinSize <- as.numeric(ctlp_df$WinSize)
ctlp_df$Chrom <- as.numeric(ctlp_df$Chrom)
ctlp_df$Start <- as.numeric(ctlp_df$Start)
ctlp_df$End <- as.numeric(ctlp_df$End)
ctlp_df$ExpSwitchNo <- as.numeric(ctlp_df$ExpSwitchNo)
ctlp_df$SwitchNo <- as.numeric(ctlp_df$SwitchNo)
ctlp_df$LR <- as.numeric(ctlp_df$LR)
ctlp_df$log_likelihoood <- log10(ctlp_df$LR)
saveRDS(ctlp_df, "output/TCGA_selected_CTLP_results.RDS")

# the default thresholds are: SwitchNo >= 20 and log10(LR) >= 8.
ct_df <- ctlp_df[ctlp_df$SwitchNo >= 20 & ctlp_df$log_likelihoood >= 8, ]

metadata_df$CT <- ifelse(metadata_df$patient %in% ct_df$patient_barcode, "CT+", "CT-")


saveRDS(metadata_df, "data/selected_data/meta.RDS")
