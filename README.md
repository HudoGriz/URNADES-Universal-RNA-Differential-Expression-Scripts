# URNADES
### Universal RNA Differential Expression Scripts

This R scripts perform differential expression analysis (DE) on RNA-seq data. They can be used with either STAR-featureCounts or SALMON data.

For example `HTML` reports please see the `example` directory.

### Usage

To run, you will need to provide the following information:

    The path to the sample information file (a .csv with at least sampleID and Treatment columns).
    The formula for the DE analysis (batch example: "~ Treatment", or paired sample: "~ PatientID + Treatment").
    The name of the condition of interest (example: "Treatment").
    The path to the STAR-featureCounts or SALMON data.
    The output path.
    The minimum count (optional, default = 10).
    The FDR value (optional, default = 0.05).
    The minimal Log2 fold change.

The STAR-featureCounts should be combined into a `counts.tsv` file (column per sample, row per gene).
Path to the SALMON files can be given directly -- the parent directory, which contains `sampleID/quant.sf`.
Make sure the `sampleID` matches the one in `sample_info.csv`.

To ensure the presence of all required R libraries simply crete a conda env from the `environment.yaml`:
```
conda env create -f environment.yaml
-or-
mamba env create -f environment.yaml

conda activate urnades
```
You can also install R packages locally with the help of `renv`. Make sure that the `renv.lock` file is present.

```
# Open R terminal in project directory
# Make sure the R library renv is installed locally

renv::restore()
```

For example, to run the script with STAR-featureCounts & SALMON data, you would use the following command:

```
Rscript URNADES.R \
--sampleInfoFilePath sample_info.csv \
--STAR-featureCountsdata $STAR-featureCounts_data_path \
--SALMONdata $SALMON_data_path \
--output $output \
--formula_input "~ PatientID + Treatment" \
--conditionName "Treatment" \
```

The `sample_info.csv` should look like this:
```
SampleID,PatientID,Treatment
C1_S3,Patient_1,Pre
D1_S4,Patient_1,Semaglutid
G1_S7,Patient_2,Pre
H1_S8,Patient_2,Semaglutid
```

### Output

The scripts will create two HTML reports in the output directory:

* `STAR-featureCounts_report.html`: This report contains the results of the DE analysis for STAR-featureCounts data.
* `SALMON_report.html`: This report contains the results of the DE analysis for SALMON data.

Additionally, a subfolder containing analysis result data tables is added.

### License

This scripts are licensed under the MIT License.
