suppressPackageStartupMessages({
	library(dplyr); library(ggvenn)
})

venn_and_common <- function(res_deseq2, res_edger, fdr, log2FCT = NULL, annotated, out_csv) {
	# sets based on adj p
	sel_deseq <- res_deseq2$padj < fdr & !is.na(res_deseq2$padj)
	sel_edger <- res_edger$FDR	 < fdr & !is.na(res_edger$FDR)

	if (!is.null(log2FCT)) {
		sel_deseq <- sel_deseq & (abs(res_deseq2$log2FoldChange) >= log2FCT)
		# edgeR: use logFC threshold
		sel_edger <- sel_edger & (abs(res_edger$logFC) >= log2FCT)
	}

	genes_DESeq2 <- res_deseq2$ID[sel_deseq]
	genes_edgeR	<- res_edger$ID[sel_edger]

	D <- list(DESeq2 = genes_DESeq2, edgeR = genes_edgeR)
	venn_plot <- ggvenn(D, show_percentage = TRUE)

	common <- intersect(genes_DESeq2, genes_edgeR)

	# build combined table
	r2 <- res_deseq2
	re <- res_edger
	rownames(r2) <- r2$ID; rownames(re) <- re$ID
	common_df <- dplyr::bind_cols(
		ID						= common,
		log2FoldChange_DESeq2 	= r2[common, "log2FoldChange", drop = TRUE],
		log2FoldChange_edgeR	= re[common, "logFC", drop = TRUE],
		pvalue_DESeq2			= r2[common, "pvalue", drop = TRUE],
		pvalue_edgeR			= re[common, "PValue", drop = TRUE],
		padj_DESeq2				= r2[common, "padj", drop = TRUE],
		padj_edgeR				= re[common, "FDR", drop = TRUE]
	)
	ann <- annotated[annotated[[names(annotated)[1]]] %in% common, , drop = FALSE]
	names(ann)[names(ann) == names(annotated)[1]] <- "ID"
	common_df <- dplyr::left_join(as.data.frame(common_df), ann, by = "ID")

	# write
	write.table(common_df, file = out_csv, sep = "\t", quote = FALSE, row.names = FALSE)

	list(venn_plot = venn_plot, common_df = common_df)
}
