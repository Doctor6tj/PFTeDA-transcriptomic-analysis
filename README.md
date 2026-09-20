# PFTeDA transcriptomic analysis

Analysis code and machine-readable results for the study of PFTeDA responses across three independent experimental study families. Version: 1.2.0.

## Contents

- `code/`: statistical models, evidence classification, dose-sensitivity synthesis and supplementary-table export.
- `data/`: sample annotations and a manifest of public input URLs, file sizes and SHA-256 checksums.
- `environment/`: software requirements and the package versions used for the reported results.
- `reference/`: the complete pathway-result table and the manuscript's Tables S1–S14, used to check a reproduced analysis.
- `download_inputs.py`, `run_analysis.py`, `verify_results.py`: download, execution and result-comparison entry points.

The source count matrices and gene-set collection are retrieved from their public sources. Model objects, gene-level results, diagnostic plots and execution logs are generated during an analysis. The archive focuses on statistical reproducibility; journal document formatting and final figure-layout tooling are outside its scope.

## Study design

| Study family | GEO accessions | Primary contexts |
| --- | --- | --- |
| 1 | GSE336490 | HIEC-6 and HaCaT |
| 2 | GSE244107 and GSE244109 | iPSC-derived hepatic and cardiac models |
| 3 | GSE145239 | Liver spheroids at 24 h and 240 h |

The primary analysis includes 211 profiles across these contexts. Each family contributes at most one directional vote per program; the two accessions in Family 2 are one study family. BEAS-2B is a supplementary sensitivity context with unresolved validity of its reported 6160 µM exposure, and HL7702 is an exploratory context that does not represent normal hepatocytes.

Family 3 uses an equal-weight summary of the five doses shared by both durations: 0.06, 0.67, 3.35, 6.7 and 17 µM, with the matching solvent controls. The 33 µM exposure occurs only at 24 h and is excluded from the primary summary contrast. Its samples remain in model estimation. Separate analyses assess all available doses and the lowest four shared doses. The all-available-dose sensitivity uses six doses at 24 h and five at 240 h, with a lowest-four-dose duration/batch interaction. The primary duration/batch interaction uses all five shared doses.

The selected Family 3 matrices contain 40 PFTeDA and 89 vehicle profiles. Companion accession GSE144775 is outside the analysis. Statistical thresholds, gene-set coverage criteria and the random seed are in `analysis_config.json`.

## Environment

Use Python 3.12 and R 4.5.3. Run the commands below from the repository directory, using a dedicated Python environment if desired:

```text
python -m pip install -r environment/requirements.txt
```

In R 4.5.3, install the required packages. The Bioconductor packages are from release 3.22:

```r
install.packages(c("BiocManager", "remotes"), repos = "https://cloud.r-project.org")
BiocManager::install(c("edgeR", "limma", "singscore"),
                    version = "3.22", ask = FALSE, update = FALSE)
required <- read.csv("environment/R_packages.csv")
cran <- required[!required$package %in% c("edgeR", "limma", "singscore"), ]
for (i in seq_len(nrow(cran))) {
  remotes::install_version(cran$package[i], version = cran$version[i],
                          repos = "https://cloud.r-project.org", upgrade = "never")
}
```

Some platforms need a compiler toolchain to install archived source packages. Package-management details are documented by [Bioconductor](https://bioconductor.org/install/). `environment/R_packages.csv` records the seven directly required versions; `environment/R_session_packages.csv` records other packages loaded during analysis. Check the direct dependencies before proceeding:

```text
Rscript --vanilla environment/check_environment.R
```

## Reproduce the analysis

```text
python download_inputs.py
python run_analysis.py
```

If `Rscript` is not on the command search path, supply its executable location using the `--rscript` option. Choose a different output directory with `--output results_repeat`. A nonempty output directory is rejected to prevent accidental reuse of old results.

The downloader retrieves eight files: the selected count matrices and sample maps from NCBI GEO, and the Hallmark 2020 gene-set collection from Enrichr. Every downloaded file must match its recorded checksum. A changed or unavailable source requires investigation; the program will stop instead of substituting a different input. Public input locations and identifiers are retained for reproducibility.

The main command fits the models, synthesizes the three dose-analysis definitions, exports Tables S1–S14, and compares the results with `reference/`. It produces:

- `outputs/pathway_results/ALL_PATHWAY_RESULTS.csv`: 3,050 pathway tests from 11 models, covering primary, sensitivity and exploratory contexts.
- `outputs/primary_common_five/`, `outputs/sensitivity_all_available/` and `outputs/sensitivity_lowest_four/`: evidence summaries for the three dose definitions.
- `outputs/tables/`: Tables S1–S14. S1 is the documented study-design table; S2–S14 are generated from the fitted results and sample annotations.
- `outputs/result_verification.json`: numerical and categorical comparison results. Numerical tolerances are relative 1e-8 and absolute 1e-10; categorical values must agree.

Model objects and gene-level tables are retained in `outputs/intermediate/` and `outputs/gene_results/`. These generated working files are excluded by `.gitignore`.

To repeat only the evidence synthesis and table export after a completed statistical run:

```text
python run_analysis.py --synthesis-only --output outputs
```

The reported primary results comprise four R2/R3 programs, including two R3 programs; 23 supported PFTeDA–PFOA comparison rows representing 20 unique programs; eight supported context-heterogeneity records; and 11 model-restricted family-program records. Epithelial–mesenchymal transition is R3 in the primary and all-available-dose analyses and R2 in the lowest-four-dose analysis.

## Data sources

The analyzed transcriptomic data are available from NCBI GEO: [GSE336490](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE336490), [GSE244107](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE244107), [GSE244109](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE244109) and [GSE145239](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145239). The Hallmark gene-set input is the `MSigDB_Hallmark_2020` collection distributed through [Enrichr](https://maayanlab.cloud/Enrichr/). The reader normalizes the source display label `Pperoxisome` to `Peroxisome`; gene membership is unchanged. Exact download URLs and checksums are in `data/input_manifest.csv`.
