suppressPackageStartupMessages({
    library(optparse)
    library(knitr)
    library(rmarkdown)
    library(msigdbr)
    # library(biomaRt)
})

# Script location
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)
print(paste("Script location:", script.basename))

# Parameters
print("Setting Parameters")

# Define the option parser
option_list <- list(
    make_option(c("--sampleInfoFilePath"), type = "character", default = NULL, help = "Path to sample info file"),
    make_option(c("--featureCounts"), type = "character", default = NULL,
        help = "Path to STAR combined tsv [$featureCounts_data_path.tsv]"),
    make_option(c("--SALMONdata"), type = "character", default = NULL,
        help = "Path to SALMON data-folder [$SALMON_data_path/$SAMPLE/quant.sf]"),
    # make_option(c("--scRNApseudobulk"), type = "character", default = NULL,
    #     help = "Path to scRNA pseudobulk counts [$scRNA_pseudobulk_data_path.tsv]"),
    make_option(c("--conditionName"), type = "character", default = NULL, help = "Condition name"),
    make_option(c("--output"), type = "character", default = NULL, help = "Output path, use full path to a directory"),

    make_option(c("--fdr"), type = "double", default = 0.05, help = "FDR value"),
    make_option(c("--log2FCT"), type = "double", default = 2, help = "log2 fold change threshold"),
    make_option(c("--min_count"), type = "integer", default = 10, help = "Minimum count"),
    make_option(c("--formula_input"), type = "character", default = NULL, help = "Formula input"),

    make_option(c("--annotatedPath"), type = "character", default = "Data/annotated.csv",
        help = "Path to annotation csv [gene_id, gene_symbol, description], gtf or gff"),
    make_option(c("--t2gPath"), type = "character", default = "Data/genes.filtered.t2g",
        help = "Path to genes.filtered.t2t"),
    make_option(c("--gene_name"), type = "character", default = "gene_id",
        help = "How are the gnes named in the counts file, gene_id, or gene_symbol"),
    make_option(c("--set_colors"), type = "character", default = "random",
        help = "1 to sample number of Hex Color Codes seperated by comma, e.g. #FF0000,#00FF00,#0000FF or random"),

    # ---- New enrichment options ----
    make_option(c("--enrichment_sources"), type = "character", default = "",
        help = "Comma-separated sources for enrichment: DESeq2, edgeR, common_sig_genes, common_sig_DE_genes"),
    make_option(c("--species"), type = "character", default = "Homo sapiens",
        help = "Species name for msigdbr (e.g. 'Homo sapiens', 'Mus musculus'). Required if enrichment is specified"),
    make_option(c("--gs_collection"), type = "character", default = "H",
        help = "MSigDB gene set collection (e.g. H, C1-C8). If not provided, use H")
)

# Parse the command line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if required arguments are empty
if (any(sapply(opt[c("sampleInfoFilePath", "conditionName", "output")], is.null))) {
    stop("Missing required arguments. Please provide values for --sampleInfoFilePath, --conditionName, --output")
}

if (is.null(opt["featureCounts"]) & is.null(opt["SALMONdata"])) {
    stop("Missing required arguments. Please provide --featureCounts or --SALMONdata")
}

# Assign the parsed values to variables
fdr <- opt$fdr
min_count <- opt$min_count
log2FCT <- opt$log2FCT

annotatedPath <- opt$annotatedPath
t2gPath <- opt$t2gPath
gene_name <- opt$gene_name

sampleInfoFilePath <- opt$sampleInfoFilePath
featureCounts <- opt$featureCounts
SALMONdata <- opt$SALMONdata
scRNApseudobulk <- opt$scRNApseudobulk
output <- opt$output

conditionName <- opt$conditionName
formula_input <- opt$formula_input

set_colors <- opt$set_colors

enrichment_species <- opt$species
gs_collection <- opt$gs_collection

# Enrichment sources
# Split enrichment sources into vector
enrichment_sources <- c("DESeq2", "edgeR", "common_sig_genes", "common_sig_DE_genes")
enrichment_sources <- sapply(enrichment_sources, function(x) {
    if (!is.null(opt$enrichment_sources)) {
        grepl(x, opt$enrichment_sources)
    } else {
        FALSE
    }
})
enrichment_sources <- names(enrichment_sources[enrichment_sources])

# If enrichment requested but species missing then stop
if (length(enrichment_sources) > 0) {
    if (is.null(enrichment_species)) {
        stop("Error: --species must be specified when running enrichment analysis.")
    }

    # print("Loading MSigDB collections")
    print(paste("Loading MSigDB collections for species:", enrichment_species))
    print(paste("Using MSigDB collection:", gs_collection))
    MSigDB <- msigdbr::msigdbr(species = enrichment_species, collection = gs_collection)
}

# Check if output directory exists
if (!dir.exists(output)) {
    dir.create(output, recursive = TRUE)
    print("Output directory created.")
} else {
    print("Output directory already exists.")
}
output <- normalizePath(output)

print("Getting annotation database")
# Check file extensions
if (grepl(".gtf$", annotatedPath)) {
    print("Annotation database: GTF")
    gtf <- rtracklayer::import(annotatedPath)

    annotated <- data.frame(
        gene_id = gtf$gene_id,
        gene_symbol = gtf$gene_name,
        description = gtf$description
        )
}

if (grepl(".gff$", annotatedPath)) {
    print("Annotation database: GFF")
    gff <- rtracklayer::import(annotatedPath)

    annotated <- data.frame(
        gene_id = paste0("gene-", gff$gene),
        gene_symbol = gff$gene,
        description = gff$description
        )
}

# BiomaRt anotation database
# Check if file exists
if (grepl(".csv$", annotatedPath)) {
    print("Annotation database: CSV")
    if (!file.exists(annotatedPath)) {
        print("Downloading annotation database")
        # Create annotated.csv
        mart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

        annotated <- biomaRt::getBM(
            attributes = c("ensembl_gene_id", "hgnc_symbol", "description"),
            # filters = "ensembl_gene_id",
            # values = results_STAR_edgeR$ID,
            mart = mart
            )

        # Save annotated as .csv
        names(annotated) <- c("gene_id", "gene_symbol", "description")
        write.csv(annotated, file = annotatedPath, row.names = FALSE)
    } else {
        annotated <- read.csv(annotatedPath, header = TRUE, stringsAsFactors = FALSE)
    }
}


# Experiment data
print("Getting experiment data")
coldata <- read.table(sampleInfoFilePath, header = TRUE, stringsAsFactors = FALSE, sep = ",")
col_names <- colnames(coldata)

# Setting design formula
print("Setting design formula")
if (is.null(formula_input)) {
    formula <- as.formula(paste("~", paste(col_names[-1], collapse = " + ")))
} else {
    formula <- as.formula(formula_input)
}

# Condition value column number
conditionN <- which(col_names == conditionName)

# Color & size
n <- length(unique(coldata[, 1])) + length(unique(coldata[, conditionN]))
if (set_colors == "random") {
    col_vector <- grDevices::colors()[grep("gr(a|e)y", grDevices::colors(), invert = T)]
    sub_col_vector <- sample(size = n, col_vector)
} else {
    sub_col_vector <- unlist(strsplit(set_colors, ","))
    sub_col_vector <- sub_col_vector[1:n]
}


print("Starting pipeline")
data_origins <- c()

# Run STAR analysis if path is given
if (!is.null(featureCounts)) {
    print("-- Running featureCounts analysis")
    rmarkdown::render(
        paste0(
            script.basename,
            "/Rmd/edgeR_DESeq2_report.Rmd"
            ),
        output_format = "html_document",
        output_file = paste0(output, "/featureCounts_report.html"),
        clean = TRUE,
        envir = new.env(),
        params = list(
            data_origin = "featureCounts",
            countFilePath = featureCounts,
            output = output,
            coldata = coldata,
            conditionN = conditionN,
            annotated = annotated,
            min_count = min_count,
            fdr = fdr,
            log2FCT = log2FCT,
            formula = formula,
            sub_col_vector = sub_col_vector,
            gene_name = gene_name
            )
        )

    data_origins <- c(data_origins, "featureCounts")
}


# Run SALMON analysis if path is given
if (!is.null(SALMONdata)) {
    print("-- Running SALMON analysis")
    rmarkdown::render(
        paste0(
            script.basename,
            "/Rmd/edgeR_DESeq2_report.Rmd"
            ),
        output_format = "html_document",
        output_file = paste0(output, "/SALMON_report.html"),
        clean = TRUE,
        envir = new.env(),
        params = list(
            data_origin = "SALMON",
            countFilePath = SALMONdata,
            output = output,
            coldata = coldata,
            conditionN = conditionN,
            annotated = annotated,
            min_count = min_count,
            fdr = fdr,
            log2FCT = log2FCT,
            formula = formula,
            t2gPath = t2gPath,
            sub_col_vector = sub_col_vector,
            gene_name = gene_name
            )
        )
    data_origins <- c(data_origins, "SALMON")
}

# Run enrichment analysis
if (length(enrichment_sources) > 0) {
    for (data_origin in data_origins) {
        for (enrichment_source in enrichment_sources) {
            print(paste("-- Running enrichment analysis for: ", data_origin, " --> ", enrichment_source))
            rmarkdown::render(
                paste0(
                    script.basename,
                    "/Rmd/enrichment_report.Rmd"
                    ),
                output_format = "html_document",
                output_file = paste0(output, "/", data_origin, "_", enrichment_source, "_enrichment_report.html"),
                clean = TRUE,
                envir = new.env(),
                params = list(
                    enrichment_source = enrichment_source,
                    data_origin = data_origin,
                    expression_results = paste0(output, "/", data_origin, "/", enrichment_source, ".csv"),
                    r_output = output,
                    MSigDB = MSigDB,
                    gs_collection = gs_collection
                )
            )
        }
    }
}

# TODO: finish implementation
# # Run scRNA pseudobulk analysis if path is given
# if (!is.null(scRNApseudobulk)) {
#     print("-- Running scRNA pseudobulk analysis")
#     print(scRNApseudobulk)
#     rmarkdown::render(
#         paste0(
#             script.basename,
#             "/Rmd/edgeR_DESeq2_report.Rmd"
#             ),
#         output_format = "html_document",
#         output_file = paste0(output, "/scRNA_pseudobulk_report.html"),
#         clean = TRUE,
#         envir = new.env(),
#         params = list(
#             data_origin = "scRNA_pseudobulk",
#             countFilePath = scRNApseudobulk,
#             output = output,
#             coldata = coldata,
#             conditionN = conditionN,
#             annotated = annotated,
#             min_count = min_count,
#             fdr = fdr,
#             log2FCT = log2FCT,
#             formula = formula,
#             sub_col_vector = sub_col_vector,
#             gene_name = gene_name
#             )
#         )
# }
