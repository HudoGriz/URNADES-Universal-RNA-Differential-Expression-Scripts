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
    make_option(c("--STARdata"), type = "character", default = NULL, help = "Path to STAR data"),
    make_option(c("--SALMONdata"), type = "character", default = NULL, help = "Path to SALMON data"),
    make_option(c("--conditionName"), type = "character", default = NULL, help = "Condition name"),
    make_option(c("--output"), type = "character", default = NULL, help = "Output path, use full path to a directory"),

    make_option(c("--fdr"), type = "double", default = 0.05, help = "FDR value"),
    make_option(c("--log2FCT"), type = "double", default = 2, help = "log2 fold change threshold"),
    make_option(c("--min_count"), type = "integer", default = 10, help = "Minimum count"),
    make_option(c("--formula_input"), type = "character", default = NULL, help = "Formula input"),

    make_option(c("--annotatedPath"), type = "character", default = "Data/annotated.csv", help = "Path to annotated.csv"),
    make_option(c("--t2gPath"), type = "character", default = "Data/genes.filtered.t2g", help = "Path to genes.filtered.t2t")
)

# Parse the command line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if required arguments are empty
if (any(sapply(opt[c("sampleInfoFilePath", "conditionName", "output")], is.null))) {
    stop("Missing required arguments. Please provide values for --sampleInfoFilePath, --conditionName, --output")
}

if (is.null(opt["STARdata"]) & is.null(opt["SALMONdata"])) {
    stop("Missing required arguments. Please provide --STARdata or --SALMONdata")
}

# Assign the parsed values to variables
fdr <- opt$fdr
min_count <- opt$min_count
log2FCT <- opt$log2FCT

annotatedPath <- opt$annotatedPath
t2gPath <- opt$t2gPath

sampleInfoFilePath <- opt$sampleInfoFilePath
STARdata <- opt$STARdata
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

# BiomaRt anotation database
# Check if file exists
print("Getting annotation database")
if (!file.exists(annotatedPath)) {
    # Create annotated.csv
    mart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

    annotated <- biomaRt::getBM(
        attributes = c("ensembl_gene_id", "hgnc_symbol", "description"),
        # filters = "ensembl_gene_id",
        # values = results_STAR_edgeR$ID,
        mart = mart
        )

    # Save annotated as .csv
    write.csv(annotated, file = annotatedPath, row.names = FALSE)
} else {
    annotated <- read.csv(annotatedPath, header = TRUE, stringsAsFactors = FALSE)
}


# Experiment data
print("Getting experiment data")
coldata <- read.table(sampleInfoFilePath, header = TRUE, stringsAsFactors = FALSE, sep = ",")
col_names <- colnames(coldata)

# Model design
print("Setting model design")
if (is.null(formula_input)) {
    formula <- as.formula(paste("~", paste(col_names[-1], collapse = " + ")))
} else {
    formula <- as.formula(formula_input)
}

design <- model.matrix(formula, coldata)
rownames(design) <- coldata[[1]]

# Condition value column number
conditionN <- which(col_names == conditionName)


print("Running DE analysis")

# DEBUGING: Save the R environment
# save.image(file = "RNADE.RData")
# base::load("RNADE.RData")

# Run STAR analysis if path is given
if (!is.null(STARdata)) {
    print("-- Running STAR-featureCounts analysis")
    rmarkdown::render(
        "Rmd/STAR-featureCounts_report.Rmd",
        output_format = "html_document",
        output_file = paste0(output, "/STAR-featureCounts_report.html"),
        clean = TRUE,
        params = list(
            countFilePath_STAR = STARdata,
            output = output,
            coldata = coldata,
            conditionN = conditionN,
            annotated = annotated,
            min_count = min_count,
            fdr = fdr,
            log2FCT = log2FCT,
            formula = formula
            )
        )
}


# Run SALMON analysis if path is given
if (!is.null(SALMONdata)) {
    print("-- Running SALMON analysis")
    rmarkdown::render(
        "Rmd/SALMON_report.Rmd",
        output_format = "html_document",
        output_file = paste0(output, "/SALMON_report.html"),
        clean = TRUE,
        params = list(
            countFilePath_SALMON = SALMONdata,
            output = output,
            coldata = coldata,
            conditionN = conditionN,
            annotated = annotated,
            min_count = min_count,
            fdr = fdr,
            log2FCT = log2FCT,
            formula = formula,
            t2gPath = t2gPath
            )
        )
}
