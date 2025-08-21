suppressPackageStartupMessages({
	library(tximport); library(edgeR); library(DESeq2)
})

sanitize_colnames <- function(df) {
	colnames(df) <- gsub("\\.", "-", colnames(df))
	df
}
	
ensure_unique_rownames <- function(df) {
	rn <- rownames(df)
	rn2 <- sub("\\.[0-9]*$", "", rn)
	dup <- duplicated(rn2)
	if (any(dup)) {
		cat("Duplicated feature names removed: ", sum(dup), "\n")
		df <- df[!dup, , drop = FALSE]
		rownames(df) <- rn2[!dup]
	} else {
		rownames(df) <- rn2
	}
	df
}

align_counts_and_metadata <- function(countdata, coldata) {
	# countdata <- sanitize_colnames(countdata)
	# Report the samples in countdata and coldata
	cat(paste("Countdata samples:", paste(colnames(countdata), collapse = "\n"), sep = "\n"))
    cat("\n")
    cat(paste("Extra samples in samplesheet:", 
    paste(setdiff(coldata[[1]], colnames(countdata)), collapse = "\n"), sep = "\n"))
    cat("\n")
    cat(paste("Extra samples in countdata:",
    paste(setdiff(colnames(countdata), coldata[[1]]), collapse = "\n"), sep = "\n"))

	# keep only samples present in both
	keep_samples <- intersect(colnames(countdata), coldata[[1]])
	coldata <- coldata[coldata[[1]] %in% keep_samples, , drop = FALSE]
	countdata <- countdata[, coldata[[1]], drop = FALSE]

	list(countdata = countdata, coldata = coldata)
}

build_design_from_formula <- function(formula, coldata) {
	coldata_transformed <- coldata
	coldata_transformed[] <- lapply(coldata_transformed, function(x) {
    	if (is.character(x)) factor(x, levels = unique(x)) else x
	})

	design <- model.matrix(as.formula(formula), data = coldata_transformed)
	rownames(design) <- coldata[[1]]
	design
}

# load_featureCounts <- function(path) {
# 	# expects a TSV like featureCounts combined; gene ids in first column
# 	df <- read.table(path, header = TRUE, sep = "\t", check.names = TRUE, stringsAsFactors = FALSE)
# 	rownames(df) <- df[[1]]
# 	df <- df[, -1, drop = FALSE]
# 	df
# }

load_salmon <- function(root_path, t2g_path) {
	# expects salmon quant.sf files in subdirectories
	salmon_files <- list.files(root_path, recursive = TRUE, pattern = "quant.sf", full.names = TRUE)
	names(salmon_files) <- basename(dirname(salmon_files))
	tx2gene <- read.delim(t2g_path, header = FALSE)[, 1:2]
	tx <- tximport(salmon_files, type = "salmon", tx2gene = tx2gene)
	as.data.frame(tx$counts)
}

load_counts_and_align <- function(data_origin, count_path, coldata, t2g_path, gene_name, outdir) {
	# browser()
	if (identical(data_origin, "SALMON")) {
		countdata <- load_salmon(count_path, t2g_path)
	} else {
		countdata <- read.table(count_path, header = TRUE, check.names = TRUE, stringsAsFactors = FALSE, sep = ifelse(grepl("\\.tsv$", count_path, ignore.case=TRUE), "\t", " "))
		if (!is.null(rownames(countdata)) && !all(rownames(countdata) %in% countdata[[1]])) {
			# if first column is not rownames, try to set
			if (colnames(countdata)[1] != "") rownames(countdata) <- countdata[[1]]
			if (colnames(countdata)[1] != "") countdata <- countdata[, -1, drop = FALSE]
		}
	}

	countdata <- sanitize_colnames(countdata)
	countdata <- ensure_unique_rownames(countdata)

	aligned <- align_counts_and_metadata(countdata, coldata)
	countdata <- aligned$countdata
	coldata	 <- aligned$coldata


	design <- build_design_from_formula(~ . + 0 + ., coldata) # will be overwritten by passed formula later if needed

	# but we want the user-specified formula; rebuild here:
	# The caller passes `formula` from params; we just return placeholder here.
	list(countdata = countdata, coldata = coldata, design = build_design_from_formula(formula, coldata))
}
