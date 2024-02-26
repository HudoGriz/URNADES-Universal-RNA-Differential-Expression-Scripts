suppressPackageStartupMessages({
    library(optparse)
    library(knitr)
    library(rmarkdown)
    # library(biomaRt)
})


# Parameters
print("Setting Parameters")

# Define the option parser
option_list <- list(
    make_option(c("--sampleInfoFilePath"), type = "character", default = NULL, help = "Path to sample info file"),
    make_option(c("--featureCounts"), type = "character", default = NULL,
        help = "Path to STAR combined tsv [$featureCounts_data_path]"),
    make_option(c("--SALMONdata"), type = "character", default = NULL,
        help = "Path to SALMON data-folder [$SALMON_data_path/$SAMPLE/quant.sf]"),
    make_option(c("--conditionName"), type = "character", default = NULL, help = "Condition name"),
    make_option(c("--output"), type = "character", default = NULL, help = "Output path, use full path to a directory"),

    make_option(c("--fdr"), type = "double", default = 0.05, help = "FDR value"),
    make_option(c("--log2FCT"), type = "double", default = 2, help = "log2 fold change threshold"),
    make_option(c("--min_count"), type = "integer", default = 10, help = "Minimum count"),
    make_option(c("--formula_input"), type = "character", default = NULL, help = "Formula input"),

    make_option(c("--annotatedPath"), type = "character", default = "Data/annotated.csv",
        help = "Path to annotation csv [gene_id, gene_symbol, description], gtf or gff"),
    make_option(c("--t2gPath"), type = "character", default = "Data/genes.filtered.t2g",
        help = "Path to genes.filtered.t2t")
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

sampleInfoFilePath <- opt$sampleInfoFilePath
featureCounts <- opt$featureCounts
SALMONdata <- opt$SALMONdata
output <- opt$output

conditionName <- opt$conditionName
formula_input <- opt$formula_input

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
col_vector <- grDevices::colors()[grep("gr(a|e)y", grDevices::colors(), invert = T)]
n <- length(unique(coldata[, 1])) + length(unique(coldata[, conditionN]))
sub_col_vector <- sample(size = n, col_vector)


print("Running DE analysis")


# DEBUGING: Save the R environment
# save.image(file = "RNADE.RData")
# base::load("RNADE.RData")


# Run STAR analysis if path is given
if (!is.null(featureCounts)) {
    print("-- Running featureCounts analysis")
    rmarkdown::render(
        "Rmd/edgeR_DEseq2_report.Rmd",
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
            sub_col_vector = sub_col_vector
            )
        )
}


# Run SALMON analysis if path is given
if (!is.null(SALMONdata)) {
    print("-- Running SALMON analysis")
    rmarkdown::render(
        "Rmd/edgeR_DEseq2_report.Rmd",
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
            sub_col_vector = sub_col_vector
            )
        )
}
