suppressPackageStartupMessages({
	library(ggplot2); library(ggvenn)
})

pretty_n <- function(x) format(x, big.mark = ",", scientific = FALSE)


# ---------------------------
# Plot log-CPM densities
# ---------------------------
plot_lcpm_density <- function(y_raw, y_filtered, L = NULL, M = NULL, main_titles = c("A. Raw data", "B. Filtered data")) {
	# Compute mean and median library size if not provided
	if (is.null(L)) L <- mean(y_raw$samples$lib.size) * 1e-6
	if (is.null(M)) M <- median(y_raw$samples$lib.size) * 1e-6

	# log-CPM cutoff
	lcpm.cutoff <- log2(10 / M + 2 / L)
	nsamples <- ncol(y_raw)
	col <- RColorBrewer::brewer.pal(nsamples, "Paired")

	par(mfrow = c(1, 2))

	# Function to plot density for one dataset
	plot_density <- function(y_data, title_text) {
		lcpm <- edgeR::cpm(y_data, log = TRUE)
		plot(density(lcpm[, 1]), col = col[1], lwd = 2, las = 2, main = "", xlab = "")
		title(main = title_text, xlab = "Log-cpm")
		abline(v = lcpm.cutoff, lty = 3)
		for (i in 2:nsamples) {
		den <- density(lcpm[, i])
		lines(den$x, den$y, col = col[i], lwd = 2)
		}
	}

	plot_density(y_raw, main_titles[1])
	plot_density(y_filtered, main_titles[2])

	invisible(dev.off())
}
