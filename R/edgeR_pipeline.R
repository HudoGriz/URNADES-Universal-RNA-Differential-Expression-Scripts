# edgeR_pipeline.R
suppressPackageStartupMessages({
    library(edgeR)
    library(ggplot2)
    library(plotly)
    library(RColorBrewer)
    library(EnhancedVolcano)
    library(DT)
})

# ---------------------------
# 1. Create DGEList
# ---------------------------
edger_create_dgelist <- function(countdata, design, coldata, conditionN) {
    y <- DGEList(
        counts = countdata,
        group = factor(coldata[[conditionN]]),
        design = design
    )
    all_unique <- length(rownames(y$counts)) == length(unique(rownames(y$counts)))
    print(paste("Are all genes unique in counts: ", all_unique))
    return(y)
}

# ---------------------------
# 2. Filter counts
# ---------------------------
edger_filter_counts <- function(y, design, min_count) {
    cpm <- cpm(y)
    lcpm <- cpm(y, log = TRUE)
    L <- mean(y$samples$lib.size) * 1e-6
    M <- median(y$samples$lib.size) * 1e-6
    print(paste("mean:", L, " --- median:", M))

    missing_gene <- table(rowSums(y$counts == 0))
    missing_sample_0 <- colSums(y$counts == 0)
    missing_sample_min <- colSums(y$counts < min_count)
    plot(missing_gene)

    keep <- filterByExpr(y, design, min.count = min_count)
    y_filtered <- y[keep, , keep.lib.sizes = FALSE]
    print(paste("Max:", length(keep)))
    print(paste("Keeping ", sum(keep), " genes out of ", length(keep)))
    print(paste("To keep [%]:", round((sum(keep) / length(keep)) * 100), "%"))

    plot_lcpm_density(y, y_filtered)

    # # log-CPM
    # lcpm.cutoff <- log2(10 / M + 2 / L)
    # nsamples <- ncol(y)
    # col <- brewer.pal(nsamples, "Paired")
    
    # par(mfrow = c(1, 2))
    # plot(density(lcpm[, 1]), col = col[1], lwd = 2, las = 2, main = "", xlab = "") #ylim = c(0, 0.26),
    # title(main = "A. Raw data", xlab = "Log-cpm")
    # abline(v = lcpm.cutoff, lty = 3)
    # for (i in 2:nsamples) {
    #     den <- density(lcpm[, i])
    #     lines(den$x, den$y, col = col[i], lwd = 2)
    # }

    # lcpm <- cpm(y_filtered, log = TRUE)
    # plot(density(lcpm[, 1]), col = col[1], lwd = 2, las = 2, main = "", xlab = "") #ylim = c(0, 0.26),
    # title(main = "B. Filtered data", xlab = "Log-cpm")
    # abline(v = lcpm.cutoff, lty = 3)
    # for (i in 2:nsamples) {
    #     den <- density(lcpm[, i])
    #     lines(den$x, den$y, col = col[i], lwd = 2)
    # }
    # invisible(dev.off())

    return(list(y = y_filtered, keep = keep))
}

# ---------------------------
# 3. Plot missingness
# ---------------------------
edger_plot_missingness <- function(y, min_count, coldata, conditionN, outdir) {
    missing_gene <- table(rowSums(y$counts == 0))
    missing_sample_0 <- colSums(y$counts == 0)
    missing_sample_min <- colSums(y$counts < min_count)

    missing_sample_0_df <- data.frame(sample = names(missing_sample_0), genes = missing_sample_0)
    missing_sample_min_df <- data.frame(sample = names(missing_sample_min), genes = missing_sample_min)

    p1 <- ggplot(missing_sample_0_df, aes(x = sample, y = genes)) +
        geom_bar(stat = "identity") +
        geom_hline(aes(yintercept = mean(genes), col = "mean")) +
        geom_hline(aes(yintercept = median(genes), col = "median")) +
        theme_bw()

    p2 <- ggplot(missing_sample_min_df, aes(x = sample, y = genes)) +
        geom_bar(stat = "identity") +
        geom_hline(aes(yintercept = mean(genes), col = "mean")) +
        geom_hline(aes(yintercept = median(genes), col = "median")) +
        theme_bw()

    ggplotly(p1)
    ggplotly(p2)
}

# ---------------------------
# 4. TMM normalization
# ---------------------------
edger_normalize_tmm <- function(y, design, fdr) {
    nsamples <- ncol(y)
    col <- brewer.pal(nsamples, "Paired")

    # 1. Determine max label length for margin adjustment
    max_label_len <- max(nchar(colnames(y$counts)))
    par(mfrow = c(1, 2), mar = c(1 + max_label_len * 0.5, 4, 4, 2) + 0.1)

    # 2. Boxplot of unnormalized data
    lcpm_raw <- cpm(y, log = TRUE)
    boxplot(lcpm_raw, las = 2, col = col, main = "")
    title(main = "Unnormalized data", ylab = "Log-CPM")

    # 3. TMM normalization
    y_norm <- calcNormFactors(y, method = "TMM")

    # 4. Boxplot of normalized data
    lcpm_norm <- cpm(y_norm, log = TRUE)
    boxplot(lcpm_norm, las = 2, col = col, main = "")
    title(main = "Normalized data", ylab = "Log-CPM")

    return(y_norm)
}

# ---------------------------
# 5. MDS plot
# ---------------------------
edger_plot_mds <- function(y, coldata, conditionN, outdir, sub_col_vector) {
    mds <- plotMDS(y, plot = FALSE)
    toplot <- data.frame(
        Dim1 = mds$x,
        Dim2 = mds$y,
        Group = as.factor(coldata[[conditionN]]),
        Sample = coldata[[1]]
    )

    gp <- ggplot(toplot, aes(Dim1, Dim2, color = Group)) +
        geom_hline(aes(yintercept = 0), colour = "grey") +
        geom_vline(aes(xintercept = 0), colour = "grey") +
        geom_point(size = 4) +
        xlab("Dimension 1") +
        ylab("Dimension 2") +
        ggtitle("MDS plot\n") +
        theme_bw() +
        geom_point(aes(toplot$Dim1, toplot$Dim2, color = toplot$Sample), size = 2) +
        scale_color_manual(values = sub_col_vector)

    ggplotly(gp)

}

# ---------------------------
# 6. Estimate dispersion
# ---------------------------
edger_estimate_dispersion <- function(y, design, outdir) {
    y <- estimateDisp(y, design, robust = TRUE)
    # png(file.path(outdir, "edgeR_BCV.png"), width = 1200, height = 600)
    plotBCV(y)
    # dev.off()
    return(y)
}

# ---------------------------
# 7. Fit model & test
# ---------------------------
edger_fit_and_test <- function(y, design, log2FCT, outdir) {
    fit <- glmQLFit(y, design)
    lrt <- glmLRT(fit)
    # png(file.path(outdir, "edgeR_MD.png"), width = 1200, height = 600)
    plotMD(lrt)
    abline(h = c(-log2FCT, log2FCT), col = "blue")
    # dev.off()
    return(lrt)
}

# ---------------------------
# 8. P-value histogram
# ---------------------------
edger_pvalue_histogram <- function(results, fdr, outdir) {
    hp <- ggplot(results) +
        geom_histogram(aes(PValue, fill = "p-value"), colour = "grey20", alpha = 0.5, stat = "bin") +
        geom_histogram(aes(FDR, fill = "padj"), colour = "grey20", alpha = 0.5, stat = "bin") +
        scale_fill_manual(
            name = "group", values = c("p-value" = "steelblue", "padj" = "grey20"), labels = c("a" = "p-value", "b" = "padj")
        ) +
        geom_vline(xintercept = fdr, colour = "red") +
        theme_bw()

    ggplotly(hp)
}

# ---------------------------
# 9. Annotate & save results
# ---------------------------
edger_annotate_results <- function(results, annotated, gene_name, fdr, lrt, outdir) {
    annotated_sub <- annotated[annotated[, gene_name] %in% row.names(results), ]
    results <- merge(results, annotated_sub, by.x = "ID", by.y = gene_name, all = TRUE)
    results <- results[order(results$FDR), ]
    results <- results[results$PValue <= fdr, ]

    write.table(
        results, file = file.path(outdir, "edgeR.csv"),
        quote = FALSE, sep = "\t", row.names = FALSE
    )

    counts <- data.frame(cpm(lrt))
    counts$ID <- rownames(counts)
    results_counts <- merge(results, counts, by = "ID", all = TRUE)
    write.table(
        results_counts, file = file.path(outdir, "edgeR_counts.csv"),
        quote = FALSE, sep = "\t", row.names = FALSE
    )

    return(results)
}

# ---------------------------
# 10. Volcano plot
# ---------------------------
edger_volcano_plot <- function(results, fdr, outdir) {
    # png(file.path(outdir, "edgeR_volcano.png"), width = 1200, height = 900, res = 150)
    EnhancedVolcano(
        results,
        lab = "",
        x = "logFC",
        y = "FDR",
        title = "edgeR results",
        legendPosition = "top",
        pCutoff = fdr
    )
    # dev.off()
}

# ---------------------------
# 11. Main wrapper
# ---------------------------
run_edger_pipeline <- function(
    countdata, design, coldata, conditionN, min_count,
    fdr, log2FCT, annotated, gene_name, outdir
) {

    y <- edger_create_dgelist(countdata, design, coldata, conditionN)
    y <- edger_filter_counts(y, design, min_count)$y
    edger_plot_missingness(y, min_count, coldata, conditionN, outdir)
    y <- edger_normalize_tmm(y, design, fdr)
    edger_plot_mds(y, coldata, conditionN, outdir)
    y <- edger_estimate_dispersion(y, design, outdir)
    lrt <- edger_fit_and_test(y, design, log2FCT, outdir)

    results <- topTags(lrt, adjust.method = "BH", n = nrow(lrt$table))
    results <- na.omit(as.data.frame(results))
    results$ID <- rownames(results)

    edger_pvalue_histogram(results, fdr, outdir)
    edger_annotate_results(results, annotated, gene_name, fdr, lrt, outdir)
    edger_volcano_plot(results, fdr, outdir)

    return(results)
}
