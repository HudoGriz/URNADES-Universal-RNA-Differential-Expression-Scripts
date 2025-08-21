suppressPackageStartupMessages({
	library(ggplot2); library(gridExtra); library(edgeR); library(RColorBrewer)
})

dt_opts_wide <- function() {
	list(
		columnDefs = list(list(className = "dt-center", targets = "_all")),
		scrollX = TRUE, pageLength = 10, language = list(search = "Filter:")
	)
}

quick_condition_means_and_scatter <- function(countdata, coldata, conditionN) {
	cond_col <- coldata[, conditionN]
	uniq <- unique(cond_col)

	gene_means <- sapply(uniq, function(cond) {
		rowMeans(countdata[, coldata[,1][cond_col == cond], drop = FALSE])
	})
	gene_means <- as.data.frame(gene_means)
	colnames(gene_means) <- make.names(uniq)

	combs <- combn(colnames(gene_means), 2, simplify = FALSE)
	grobs <- list()
	for (i in seq_along(combs)) {
		x <- combs[[i]][1]; y <- combs[[i]][2]
		p <- ggplot(gene_means, aes_string(x = x, y = y)) +
			geom_point(alpha = 0.5) +
			scale_x_continuous(trans = "log10") +
			scale_y_continuous(trans = "log10") +
			labs(x = x, y = y) + theme_bw()
		grobs[[i]] <- p
	}

	if (length(grobs) == 1) {
		grob_plot <- grobs[[1]]
	} else {
		grob_plot <- do.call(gridExtra::arrangeGrob, c(grobs, ncol = 2))
	}

	list(gene_means = gene_means, mean_combinations = combs, grid = grob_plot)
}

add_pairwise_log2fc <- function(gene_means, mean_combinations, log2FCT) {
	gm <- gene_means
	grobs <- list()
	for (i in seq_along(mean_combinations)) {
		x <- mean_combinations[[i]][1]; y <- mean_combinations[[i]][2]
		colname <- paste0("log2FC_", x, "_", y)
		gm[[colname]] <- log2(gm[[y]] / gm[[x]])

		# Report
		print(paste0(colname, " - Differently expressed: ", sum(abs(gm[[colname]]) > log2FCT)))
		print(paste0(colname, " - Underexpressed: ", sum(gm[[colname]] < -log2FCT)))
		print(paste0(colname, " - Overexpressed: ", sum(gm[[colname]] > log2FCT)))

		sig <- abs(gm[[colname]]) > log2FCT
		p <- ggplot(gm, aes_string(x = x, y = y)) +
			geom_point(alpha = 0.3) +
			geom_point(data = gm[sig, , drop = FALSE], aes(col = "red"), alpha = 0.9) +
			scale_x_continuous(trans = "log10") + scale_y_continuous(trans = "log10") +
			theme_bw() + ggtitle(colname)
		grobs[[i]] <- p
	}
	if (length(grobs) == 1) {
		grob_plot <- grobs[[1]]
	} else {
		grob_plot <- do.call(gridExtra::arrangeGrob, c(grobs, ncol = 2))
	}

	list(gene_means_log = gm, plot_grid = grob_plot)
}

