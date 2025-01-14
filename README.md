# cfDNA in Disease Biology
Source Code of Master's Thesis "XX" by Timo Rieger, February 2025 completed at Krauthammers Lab.

This git repository provides scripts and detials for the analyses performed in the manuscript "XX".

_All scripts are provided as is, without any warranty. This description and code is provided in addition to the methods in the original work and should simplify the reproduction of the analyses._

## Version and Dependencies

## Primary Data Processing - Snyder et al. 2016
cfDNA data was pre-processed as described in Snyder et al. 2016.

Link to their Git Repo: [ShendurLab cfDNA](https://github.com/shendurelab/cfDNA/blob/master)

Below are the names of the original scripts that were used and modified in this work, ordered chronologically. They can be found on their Git Repo in the folder "expression".

1. extractReadStartsFromBAM_Region_WPS.py
2. fft_path.R
3. convert_files.py
4. plots.R


### Consensus Coding Sequence - Gene Annotation File
Gene annotation was performed using the consensus coding sequence from Gencode, Gencode Reference 30 (GRCh38) genome. The **main annotation file** in **GFF3** format was downloaded from the official GENCODE website: [GENCODE 30v](https://www.gencodegenes.org/human/release_30.html)

The original file was modifed to suit the source code from Snyder et al. (2016). The resulting gene annotation file contains 5 columns: ENSG#, chromosome#, start, end and strand (+ or -).

`zcat gencode.v30.annotation.gff3.gz | awk -F '\t' '$3 == "gene" {split($9, a, ";"); split(a[1], b, "="); split(b[2], c, "."); gsub("chr", "", $1); print c[1], $1, $4, $5, $7}' OFS='\t' > gencode.v30.annotation.tsv`

Only active genes were used, thus inactive and long non-coding regions were filtere out. The final gene annotation file contained 31'591 genes (58'870 - 13'160 (inactive regions) - 13'759 (lnc regions)).

Furthermore, the frist 10 kb of each gene (TSS) was used for downstream analysis. This was done by modifing the gene annotation file accordingly:

`cat filtered_annotation_file.tsv | awk 'BEGIN{ FS="\t"; OFS="\t" } NR > 1 { if ($5 == "+") { print $1,$2,$3-1,$3-1+10000,$5 } else { print $1,$2,$4-1-10000,$4-1,$5 } }' > final_annotation_file.tsv`


## Gene Expression Reference Matrix - Tabula Sapiens
The Single-cell RNA sequencing data that was downloaded can be found under this link: [Tabula Sapiens Data](https://cellxgene.cziscience.com/collections/e5f58829-1a66-40b5-a624-9046778e74f5). Specifically, the RDS file called "Tabula Sapiens - All cells" was downloaded. More information and tools to visualize the data beforehand can be found under the offical Tabula Sapiens link: [Offical Tabula Sapiens Website](https://tabula-sapiens.sf.czbiohub.org/).

The data was loaded into R V(4:XX) and subset using the "10 x 3" assay (n=412'848 cells). Sperm cells were excluded and an additional metadata column called "cell_type_tissue" was created by concatenating the columns "cell_type" and "tissue_in_publication" resulting in 447 unique identifiers accross 24 biopsied organs. 



### Adding Placenta Data

## ATAC Seq Reference Matrix - 

## Predicitive Analyses - Stanley et al. 2024

## Evaluations
