# URNADES - Universal RNA Differential Expression Scripts

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R](https://img.shields.io/badge/R-4.3%2B-blue.svg)](https://www.r-project.org/)
[![Singularity](https://img.shields.io/badge/Singularity-Ready-green.svg)](https://sylabs.io/singularity/)

A comprehensive R pipeline for RNA differential expression analysis supporting multiple input formats and enrichment analysis.

## 🧬 Overview

URNADES is a flexible and automated pipeline for RNA-seq differential expression analysis that supports:

- **Multiple input formats**: featureCounts, SALMON quantification
- **Dual analysis engines**: DESeq2 and edgeR
- **Gene set enrichment analysis**: Integration with MSigDB collections
- **Automated reporting**: HTML reports with interactive visualizations
- **Containerization**: Ready-to-use Singularity container

## 🛠 Installation

### Option 1: Singularity Container (Recommended)

Build the container:
```bash
sudo singularity build urnades.sif Singularity.def
```

Pull from a registry:
```bash
singularity pull urnades.sif library://blazv/urnades/urnades
```

### Option 2: Manual Installation

Install required R packages:
```r
# Install BiocManager
install.packages("BiocManager")

# CRAN packages
install.packages(c(
  "optparse", "knitr", "rmarkdown", "tidyverse", 
  "ggplot2", "plotly", "RColorBrewer", "DT", 
  "gridExtra", "ggvenn", "calibrate", "msigdbr"
))

# Bioconductor packages
BiocManager::install(c(
  "DESeq2", "edgeR", "tximport", "EnhancedVolcano", "fgsea"
))
```

## 🚀 Quick Start

### Basic Usage

```bash
# Using Singularity container
singularity exec urnades.sif Rscript URNADES.R \
  --sampleInfoFilePath samples.csv \
  --featureCounts counts.tsv \
  --conditionName treatment \
  --output results/

# Direct R execution
Rscript URNADES.R \
  --sampleInfoFilePath samples.csv \
  --SALMONdata salmon_quants/ \
  --conditionName treatment \
  --output results/
```

### Input Files Required

1. **Sample Information File** (`samples.csv`):
```csv
sample,treatment,batch
sample1,control,1
sample2,control,1
sample3,treated,1
sample4,treated,1
```

2. **Count Data** (one of):
   - featureCounts output (`counts.tsv`)
   - SALMON quantification directory

## 📊 Usage Examples

### 1. Basic DESeq2/edgeR Analysis

```bash
singularity exec urnades.sif Rscript URNADES.R \
  --sampleInfoFilePath metadata.csv \
  --featureCounts gene_counts.tsv \
  --conditionName condition \
  --output analysis_results/ \
  --fdr 0.05 \
  --log2FCT 1.5
```

### 2. SALMON Data with Enrichment Analysis

```bash
singularity exec urnades.sif Rscript URNADES.R \
  --sampleInfoFilePath samples.csv \
  --SALMONdata salmon_output/ \
  --conditionName treatment \
  --output results/ \
  --enrichment_sources "DESeq2,edgeR,common_sig_DE_genes" \
  --species "Homo sapiens" \
  --gs_collection "H"
```

### 3. Custom Design Formula (paired samples)

```bash
singularity exec urnades.sif Rscript URNADES.R \
  --sampleInfoFilePath complex_design.csv \
  --featureCounts counts.tsv \
  --conditionName treatment \
  --formula_input "~ Patient_ID + Condition" \
  --output results/
```

## 🔧 Parameters

### Required Parameters

| Parameter | Description |
|-----------|-------------|
| `--sampleInfoFilePath` | Path to sample metadata CSV file |
| `--conditionName` | Column name for the main condition |
| `--output` | Output directory path |

### Data Input (at least one required)

| Parameter | Description |
|-----------|-------------|
| `--featureCounts` | Path to featureCounts/STAR output TSV |
| `--SALMONdata` | Path to SALMON quantification directory |

### Analysis Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--fdr` | 0.05 | False discovery rate threshold |
| `--log2FCT` | 2 | Log2 fold change threshold |
| `--min_count` | 10 | Minimum count for gene filtering |
| `--formula_input` | Auto | Custom design formula |

### Enrichment Analysis

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--enrichment_sources` | None | Sources for enrichment (DESeq2,edgeR,common_sig_DE_genes) |
| `--species` | "Homo sapiens" | Species for MSigDB |
| `--gs_collection` | "H" | MSigDB collection (H, C1-C8) |

### Annotation

If left unspecified, human references will be downloaded from **Ensembl**: `biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")`

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--annotatedPath` | Data/annotated.csv | Gene annotation file (CSV/GTF/GFF) |
| `--t2gPath` | Data/genes.filtered.t2g | Transcript-to-gene mapping |
| `--gene_name` | gene_id | Gene identifier type |

## 📈 Output Files

The pipeline generates:

### Analysis Reports
- `featureCounts_report.html` / `SALMON_report.html` - Main analysis reports
- `[source]_enrichment_report.html` - Enrichment analysis reports

### Data Files
```
output/
├── [data_origin]/
│   ├── DESeq2.csv              # DESeq2 results
│   ├── edgeR.csv               # edgeR results  
│   ├── common_sig_genes.csv    # Overlapping significant genes
│   ├── common_sig_DE_genes.csv # Overlapping DE genes
│   ├── DESeq2_counts.csv       # Normalized counts
├── *.html                      # Generated reports
```

## 🧪 Example Workflow

1. **Prepare your data**:
   ```bash
   # Sample metadata
   echo "sample,condition,batch" > samples.csv
   echo "ctrl1,control,1" >> samples.csv
   echo "ctrl2,control,2" >> samples.csv
   echo "treat1,treatment,1" >> samples.csv
   echo "treat2,treatment,2" >> samples.csv
   ```

2. **Run the analysis**:
   ```bash
   singularity exec urnades.sif Rscript URNADES.R \
     --sampleInfoFilePath samples.csv \
     --featureCounts gene_counts.tsv \
     --conditionName condition \
     --output my_analysis/ \
     --enrichment_sources "DESeq2,common_sig_DE_genes" \
     --species "Homo sapiens"
   ```

3. **View results**:
   Open `my_analysis/featureCounts_report.html` in your browser

## 🔬 Gene Set Collections

Available MSigDB collections:

| Collection | Description |
|------------|-------------|
| **H** | Hallmark gene sets |
| **C1** | Positional gene sets |
| **C2** | Curated gene sets |
| **C3** | Regulatory target gene sets |
| **C4** | Computational gene sets |
| **C5** | Ontology gene sets |
| **C6** | Oncogenic signature gene sets |
| **C7** | Immunologic signature gene sets |
| **C8** | Cell type signature gene sets |

## 🐛 Troubleshooting

### Common Issues

1. **Missing packages**: Use the Singularity container for guaranteed compatibility
2. **Memory issues**: Increase available RAM or filter low-count genes more aggressively
3. **File format errors**: Ensure CSV files use comma separators and proper headers

### Getting Help

- Check the generated HTML reports for diagnostic information
- Verify input file formats match expected structure
- Use `--help` flag to see all available options

## 🤝 Contributing

We welcome contributions! Please:

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Add tests if applicable
5. Submit a pull request

## 📝 Citation

If you use URNADES in your research, please cite:

```
URNADES: Universal RNA Differential Expression Scripts
Vrhovšek et al. (2025)
GitHub: https://github.com/HudoGriz/URNADES-Universal-RNA-Differential-Expression-Scripts
```

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙋 Support

For questions and support:

- 📧 Open an issue on GitHub
- 📖 Check the documentation in the HTML reports
- 🔍 Search existing issues for solutions

---

**Happy analyzing! 🧬📊**
