# cfDNA in Disease Biology
Source Code of the Master's Thesis **"Integrative Analysis of Cell-free DNA Fragmentation and Single-Cell Omics in Disease Biology"** by _Timo Rieger_, February 2025 completed in Krauthammer's Lab at University Hospital of Zurich.

This git repository provides scripts and detials for the analyses performed in the manuscript.

_All scripts are provided as is, without any warranty. This description and code is provided in addition to the methods in the original work and should simplify the reproduction of the analyses._

## Primary Data Processing - Snyder et al. (2016)
### Window Protection Score (WPS) - Feature Extraction
To extract WPS from cfDNA data, the BAM files were pre-processed as described in Snyder et al..

Link to their Git Repo: [ShendurLab cfDNA](https://github.com/shendurelab/cfDNA/blob/master)

Below are the names of the original scripts that were used and modified in this work, ordered chronologically. They can be found on their Git Repo in the folder "expression".

1. extractReadStartsFromBAM_Region_WPS.py
2. fft_path.R
3. convert_files.py
4. plots.R

The scripts used in this work regarding preprocessing for WPS features can be found under "_preprocessing/WPS_".

### Fragment Center Count (FCC) - Feature Extraction
To extract the FCC, the first script was modified. Afterwards, the signal was smoothed and the mean was calculated to assign ranks. The respective scripts are in this folder "_preprocessing/FCC_".

### Consensus Coding Sequence - Gene Annotation File
Gene annotation was performed using the consensus coding sequence from Gencode, Gencode Reference 30 (GRCh38) genome. The **main annotation file** in **GFF3** format was downloaded from the official GENCODE website: [GENCODE 30v](https://www.gencodegenes.org/human/release_30.html)

The original file was modifed to suit the source code from Snyder et al.. The resulting gene annotation file contains 5 columns: ENSG#, chromosome#, start, end and strand (+ or -).

`zcat gencode.v30.annotation.gff3.gz | awk -F '\t' '$3 == "gene" {split($9, a, ";"); split(a[1], b, "="); split(b[2], c, "."); gsub("chr", "", $1); print c[1], $1, $4, $5, $7}' OFS='\t' > gencode.v30.annotation.tsv`

Only active genes were used, thus inactive and long non-coding regions were filtere out. The final gene annotation file contained 31'591 genes (58'870 - 13'160 (inactive regions) - 13'759 (lnc regions)).

Furthermore, the frist 10 kb of each gene (TSS) was used for downstream analysis. This was done by modifing the gene annotation file accordingly:

`cat filtered_annotation_file.tsv | awk 'BEGIN{ FS="\t"; OFS="\t" } NR > 1 { if ($5 == "+") { print $1,$2,$3-1,$3-1+10000,$5 } else { print $1,$2,$4-1-10000,$4-1,$5 } }' > final_annotation_file.tsv`

The final gene annotation file used for the analysis is in the data folder: "_data_".

## Gene Expression Reference Matrix - Tabula Sapiens
The Single-cell RNA sequencing data that was downloaded can be found under this link: [Tabula Sapiens Data](https://cellxgene.cziscience.com/collections/e5f58829-1a66-40b5-a624-9046778e74f5). Specifically, the RDS file called "Tabula Sapiens - All cells" was downloaded. More information and tools to visualize the data beforehand can be found under the offical Tabula Sapiens link: [Offical Tabula Sapiens Website](https://tabula-sapiens.sf.czbiohub.org/).

The data was loaded into R V(4:XX) and subset using the "10 x 3" assay (n=412'848 cells). Sperm cells were excluded and an additional metadata column called "cell_type_tissue" was created by concatenating the columns "cell_type" and "tissue_in_publication" resulting in 447 unique identifiers accross 24 biopsied organs.

The script to filter the Tabula Sapiens data can be found under "_preprocessing/00_preprocessing_Tabula_Sapiens_filtered.R".

## Adding Placenta Cells
Placenta cells were added to the Tabula Sapiens reference dataset to perform the analysis on another dataset containing preeclampsia samples collected from the University Hospital of Zurich.

Two datasets containing placenta specific cells were downloaded:
 - [Placenta Dataset: ArrayExpress (E-MTAB-6701)](https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-6701)
- [Fetal Cell Atlas](https://descartes.brotmanbaty.org/bbi/human-gene-expression-during-development/)

The scripts regarding placenta data integration and preeclampsia classifier can be found on the following folder: "_preeclampsia_".

## ATAC-Seq Reference Matrix
For the comparison of different reference data types, the ATAC-Seq data from [Zhang et al. (2021)](https://pubmed.ncbi.nlm.nih.gov/34774128/) was used.
The necessary scripts are under "_preprocessing/ATAC_vs_RNA_".

## Statistical and Predictive Analysis
For the statistical and predictive analysis, some scripts from Stanley et al.'s study were used and modified.

Link to their Git Repo: [Cell type signatures in cell free DNA fragmentation profiles reveal disease biology](https://github.com/JorisVermeeschLab/cfDNA_cell_of_origin?tab=readme-ov-file#cell-type-signatures-in-cell-free-dna-fragmentation-profiles-reveal-disease-biology)

The final ranks across features for all metrics x reference dataset combinations can be found in "_data/ranks_" and are used for downstream analysis.
In the folder "_cancer_", all scripts regarding statistical and predicitve analysis for cancer can be found.

They include:
- correlation value distribution across metrics x reference dataset combinations
- unsupervised clustering UMAP
- rank distribution across compartments and tissue types
- Wilcoxon rank sum test
- LOOSVM for all metrics x reference dataset combinations
- Leave-one-out OvR SVM for FCC x RNA cell type tissue
- performance evaluation of the SVMs

The results of the performance evaluation are stored in this folder: "_data/evaluation_scores_".

The average feature weights of the primary SVM model FCC x RNA cell type tissue can be found in this folder "_data/feature_weights_".

## Version and Dependencies
The following Python (v3.12.8) and R (v4.4.1) packages were used in this thesis:

### Python Packages
- pandas 2.2.3
- matplotlib 3.10.0
- scikit-learn 1.6.1
- seaborn 0.13.2

### R Packages
- Seurat 4.4.0
- SeuratObject 5.0.2
- Matrix 1.6.5
- cmake 3.29.4
- nloptr 2.1.1
- lme4 1.1.36
- pbkrtest 0.5.3
- car 3.1.3
- rstatix 0.7.2
- devtools 2.4.5
- e1071 1.7.16
- pROC 1.18.5
- fmsb 0.7.6
- umap R
- uwot 0.2.2
- viridis 0.6.5
- pheatmap 1.0.12
