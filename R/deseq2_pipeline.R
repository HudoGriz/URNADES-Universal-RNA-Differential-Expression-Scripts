# DESeq2_pipeline.R
suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
    library(plotly)
    library(RColorBrewer)
    library(EnhancedVolcano)
    library(DT)
})

# ---------------------------
# 1. Create DESeqDataSet
# ---------------------------
deseq2_create_dds <- function(countdata, coldata, design) {
    # Cut max count data
	countdata_new <- countdata
	countdata_new[countdata_new >= .Machine$integer.max] <- .Machine$integer.max

	dds <- DESeqDataSetFromMatrix(
		countData = round(countdata_new),
		colData = coldata,
		design = design
	)

    all_unique <- length(rownames(dds)) == length(unique(rownames(dds)))
    print(paste("Are all genes unique in counts: ", all_unique))
    
	return(dds)
}

# ---------------------------
# 2. Filter counts
# ---------------------------
deseq2_filter_counts <- function(dds, min_count) {
    keep <- rowSums(counts(dds) >= min_count) >= 1
    dds_filtered <- dds[keep, ]
    print(paste("Keeping ", sum(keep), " genes out of ", length(keep)))
    print(paste("To keep [%]:", round((sum(keep)/length(keep))*100), "%"))

    # log-transformed counts for density plot
	deseq_plot_lcpm_density <- function(dds, title) {
	    lcpm <- log2(counts(dds, normalized = FALSE) + 1)

		nsamples <- ncol(dds)
		col <- brewer.pal(nsamples, "Paired")

		for (i in 1:nsamples) {
			den <- density(lcpm[, i])
			if (i == 1) plot(den, col = col[i], lwd = 2, main=title)
			else lines(den$x, den$y, col = col[i], lwd = 2)
		}
	}

	par(mfrow = c(1,2))
	deseq_plot_lcpm_density(dds, "Raw data density")
	deseq_plot_lcpm_density(dds_filtered, "Filtered data density")
	invisible(dev.off())

    return(list(dds = dds_filtered, keep = keep))
}
