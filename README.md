# Targeted transcriptomic and proteomic analyses for *Pathogenesis and natural history of the Bundibugyo species of Orthoebolavirus in nonhuman primates*

> Manuscript not yet published; under peer review.
 
## Data availability

Raw Nanostring RCC files are available via NCBI GEO (ACCESSION PENDING). Nanostring thresholded count matrices are available via NCBI GEO or in [`data.xlsx`](data.xlsx).

Estimated protein concentrations from LEGENDplex assays are in [`data.xlsx`](data.xlsx).

Code to generate figures is in [`analysis.r`](analysis.r), and output from this script are in [`analysis/`](analysis/).

## Methods and preprocessing

### Targeted transcriptomics of circulating blood (Nanostring)

The expression of 770 host mRNAs was quantified from whole blood RNA via the [Nanostring NHP Immunology v2 panel](https://nanostring.com/products/ncounter-assays-panels/ncounter-rna-assays/immunology/nhp-immunology/) ([PMID 29116224](https://pubmed.ncbi.nlm.nih.gov/29116224/); Bruker #115000276) according to the manufacturer’s instructions. Raw RCC files were loaded into [nSolver v4.0](https://nanostring.com/products/ncounter-analysis-system/ncounter-analysis-solutions/), and background thresholding was performed using the default parameters. Samples that failed the nSolver internal quality checks were removed. Thresholded count matrices were exported from nSolver and analyzed with [limma v3.66.0](https://bioconductor.org/packages/release/bioc/html/limma.html) in R v4.6.0. mRNAs with an FDR-adjusted p-value and a log2 fold change >1 or <-1 were considered significantly differentially expressed. 

### Targeted proteomics from sera (LEGENDplex)

Circulating chemokines, cytokines, and other protein markers associated with inflammation were measured via [LEGENDplex bead-based multiplex immunoassay panels](https://www.biolegend.com/en-us/immunoassays/legendplex). Gamma-irradiated serum samples were assessed in duplicate for NHP Inflammation V02 (#741491, 1:4 dilution), Human Fibrinolysis (#740761, 1:40,000 dilution), Human Vascular Inflammation (#740590, 1:1000 dilution), Human Thrombosis (#740892, 1:50 dilution), and NHP Chemokine/Cytokine (#740388, 1:4 dilution) panels according to manufacturer instructions. Assay standards were prepared in batches and aliquoted across all plates to ensure batch-to-batch consistency. All optional wash steps were incorporated into the workflow to reduce background signal. Assay samples were analyzed on an Accuri C6 Plus flow cytometer (BD Biosciences). The raw `.fcs` files from each assay plate were imported into [LEGENDplex Qognit](https://legendplex.qognit.com/user/login?next=home) cloud-based Data Analysis Software Suite, which determined analyte concentrations in experimental samples via 5-parameter logistic regression curve fitting to each assay standard curve. Analyte concentration data were exported from Qognit and analyzed in R using [limma v3.66.0](https://bioconductor.org/packages/release/bioc/html/limma.html) in R v4.6.0.

![Main figure](analysis/figure-main.pdf)
