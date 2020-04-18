cellatacUtils
================

R package for handling [cellatac
pipeline](https://github.com/cellgeni/cellatac) output

## Installation

``` r
devtools::install_github("emdann/cellatacUtils")
```

## Example usage

Load pipeline output as Seurat object

``` r
library(cellatacUtils)
results.dir <- "path/to/results"
cellatac.seu <- load_cellatac_seurat(results.dir)
```

Calculate QC metrics for peaks

``` r
## Load Ensembl based annotation from Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("EnsDb.Hsapiens.v86")
library(EnsDb.Hsapiens.v86)

## Load Signac for ENCODE blacklist annotation
library(Signac)

cellatac.seu <- compute_peakQC(cellatac.seu, 
                               cluster_col = "seurat_clusters", 
                               EnsDb.genome = EnsDb.Hsapiens.v86, 
                               blacklist_gr = blacklist_hg38)
```

Calculated QC metrics:

  - tot\_count: total count of fragments overlapping peak
  - tot\_cells: total number of cells in which peak is accessible
  - frac\_cells: fraction of cells in which peak is accessible
  - max\_frac: maximum percent of cells with accessible peak in any
    given cluster (to filter out peaks that are not observed to be
    accessible in at least n% cells of cells in at least one cluster)
  - peak\_width: width of peak in bps
  - exon: is peak overlapping exon?
  - intron: is peak overlapping intron?
  - promoter: is peak overlapping promoter?
  - annotation: annotation of overlapping genomic region (promoter,
    exon, intron, intergenic)
  - gene\_name: name of overlapping gene (downstream gene for promoter
    region, NA for intergenic peaks)
  - gene\_id: ENSEMBL id of overlapping gene (downstream gene for
    promoter region, NA for intergenic peaks)
  - tss\_distance: distance to closest Transcription Start Site in bps
  - ENCODE\_blacklist: is peak overlapping ENCODE blacklist region?
