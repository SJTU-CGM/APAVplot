## 1. Introduction
 
### 1.1 Pan-genome and PAV analysis

&emsp;&emsp;Pan-genome is the collective whole-genome sequences of a given population, revealing the diversity and functional potential within that population. The PAV(Presence/absence variation) analysis is an essential step in pan-genome analysis. The core genome contains genomic regions shared by all individuals and the distributed genome is not shared by all. The distributed genome can be further divided into genomic regions shared in most members (soft-core genome), regions shared between some members (distributed/accessory genome), and regions present in only one member (unique/private genome).


### 1.2 The functions in APAVplot

&emsp;&emsp;APAVplot is a R package designed for CPAV for the subsequent analysis and visualization of PAV profile. It is efficient to explore and visualize the complex results in PAV analysis. It provides the following modules:

* Visualization of coverage : First you need to build a COV class using `get_cov_obj()`. `cov_heatmap()` shows coverage profile in a heat map. `cov_density()` shows coverage distribution of interested regions.

* PAV statistics and analysis : First you need to build a PAV class using `get_pav_obj()`. `pav_heatmap()` shows PAV heat map.  `pav_hist()`, `pav_halfviolin()` and `pav_stackbar()` present the basic statistics. `pav_cluster()` clusters the samples based on PAV table. `pav_pca()` do the PCA analysis.

* Phenotype association and visualization : `pheno_stat` performs phenotype association calculations. `pheno_heatmap()`, `pheno_manhattan()`, `pheno_block()`, `pheno_bar()` and `pheno_violin()` display the results. 

* Drawing growth curve : `sim_plot()` is used for visualization of pan/core/private genome size estimation by simulation. The input is the output table of `CPAV sim`.

* Visualization of elements : `plot_ele_cov()`, `plot_ele_pav()` and `plot_ele_depth()` are used to check elements at coverage, PAV and depth level. 


## 2. Installation

### 2.1 Installing R/RStudio
&emsp;&emsp;In order to run APAVplot, we need the R software environment, which can be freely downloaded as follows:

* Install [R](https://www.r-project.org/)
* Install [RStudio](https://www.rstudio.com/)

### 2.2 Check or install packages

```{r eval=FALSE}
packages <- c("data.table", "ggdendro", "ggplot2", "ggrepel", "ggsignif", "ggnewscale", "patchwork", "snowfall", "circlize")
lapply(packages, function(x) {
	if(!require(x, character.only = TRUE)) {
		install.packages(x, dependencies = TRUE)
	}})
if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
BiocManager::install("ComplexHeatmap")

```

### 2.3 Install metaFunc from github.

```{r eval=FALSE}
if (!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
library(devtools)
install_github("SJTU-CGM/APAVplot", build_vignettes = TRUE)
```
