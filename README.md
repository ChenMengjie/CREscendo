## CREscendo

scATAC-seq captures fragments through Tn5 tagmentation, which allows for precise analysis of regulatory elements by examining the Tn5 cleavage sites. Unfortunately, the detailed base-resolution information is often condensed and simplified in current peak-based scATAC-seq analyses. The peak calling software, such as Cellranger, generates a broad range of peak sizes, while MACS2 can identify significantly smaller peaks. However, even these smaller peaks generally exceed the size of a single cis-regulatory element (CRE) as annotated by ENCODE. To address this issue, we introduce the CREscendo framework, which is designed to analyze differential usage of CREs by examining differences in Tn5 cleavage frequencies across cell types. This approach helps to recover regulatory signals that peak-based analyses often overlook.

Specifically, CREscendo segments each peak into designated regions such as CRE1, CRE2, and so forth, along with non-CRE areas. For each peak, a frequency vector is constructed, and a chi-square test is conducted to detect statistically significant differences in cleavage frequencies across different cell types. The resulting p-values are adjusted for multiple testing, and the significant CREs are identified. The CREscendo framework is implemented in the R package,  available on GitHub.

**CREscendo** is an implementation for differential usage analysis for peaks identified from scATACseq data.


## Installation

**CREscendo** can be installed from github directly as follows:

```r
devtools::install_github("ChenMengjie/CREscendo")
```

## Tutorial


[Tutorial 1: Analysis of 10X Genomics human PBMC Single-Cell Multiome Data types](https://htmlpreview.github.io/?https://github.com/ChenMengjie/CREscendo/blob/main/Tutorials/CREscendo_10X_multiome.html)

[Pre-processed dataset 1](https://drive.google.com/file/d/1N34NZYH1LfDs4RPE9y_FGDMc7Cx2cMVj/view?usp=share_link)

[Tutorial 2: Analysis of 10X Genomics human PBMC ATAC v2 Chromium X type](https://htmlpreview.github.io/?https://github.com/ChenMengjie/CREscendo/blob/main/Tutorials/CREscendo_10X_scATAC_Chromium_X2.html)

[Pre-processed dataset 2](https://drive.google.com/file/d/1CtnCInRLnGcnz6wTp_hXPSpP9ZuMs_93/view?usp=share_link)

[Tutorial 3: Analysis of 10X Genomics adult mouse contex ATAC v2 Chromium X type](https://htmlpreview.github.io/?https://github.com/ChenMengjie/CREscendo/blob/main/Tutorials/CREscendo_10X_adult_mouse.html)

[Pre-processed dataset 3](https://drive.google.com/open?id=17QMX7Db_8LEN_oM6Hsu4hRLPzjuur1tT&usp=drive_fs)

[Human hg38 CRE annotations](https://github.com/ChenMengjie/CREscendo/blob/main/data/GRCh38-cCREs.bed)


### Author

**Mengjie Chen** (U Chicago)

Bug report, comments or questions please send to mengjiechen@uchicago.edu.
