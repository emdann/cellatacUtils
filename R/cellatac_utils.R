########################
#### cellatac utils ####
########################

# library(tidyverse)
# library(Seurat)
# library(Signac)
# library(GenomicRanges)
# library(EnsDb.Hsapiens.v86)

#### Loading cellatac output ####

#' Load peaks x cell matrix
#'
#' @param out.dir path of cellatac output directory
#' @param as.Seurat logical indicating whether output should be a SeuratAssay object. If FALSE, a sparseMatrix is returned (default: TRUE)
#'
#' @return peak matrix (SeuratAssay object or sparseMatrix)
#'
#' @import Matrix
#' @import Seurat
#'
#' @export
load_cellatac_peakmat <- function(out.dir, as.Seurat=TRUE){
  mat_file <- paste0(out.dir, "peak_matrix/peaks_bc_matrix.mmtx.gz")
  peaks_file <- paste0(out.dir, "peak_matrix/peaks.txt")
  bc_file <- paste0(out.dir, "peak_matrix/bc.txt")

  mat <- readMM(mat_file)
  peaks <- scan(peaks_file, "")
  bc <- scan(bc_file, "")

  rownames(mat) <- peaks
  colnames(mat) <- bc

  if (as.Seurat){
    peak_mat_obj <- CreateAssayObject(mat, min.cells = 3, min.features = 1)
  } else {
    peak_mat_obj <- mat
  }
  return(peak_mat_obj)
}

#' Load cellatac output as Seurat object
#'
#' @param out.dir path of cellatac output directory
#'
#' @return SeuratObject with 2 assays (windows and peaks)
#'
#' @import Matrix
#' @import Seurat
#'
#' @export
load_cellatac_seurat <- function(out.dir){
  ## Load Seurat object (windows x cell)
  win_seurat <- readRDS(list.files(paste0(out.dir, "qc"), pattern="rds", full.names = TRUE))
  ## Save reductions separately because RenameAssays removes them (should be fixed in new commit)
  win_seurat <- RenameAssays(win_seurat, peaks='win')
  names(win_seurat@reductions) <- paste0(names(win_seurat@reductions), "_win")

  ## Load peaks x cell matrix and add to seurat obj
  peaks_mat <- load_cellatac_peakmat(out.dir)
  win_seurat@assays[["peaks"]] <- peaks_mat
  DefaultAssay(win_seurat) <- "peaks"

  return(win_seurat)
  }

#### Peaks QC ####

#' Compute QC metric for peaks
#'
#' @param win_seu Seurat object
#' @param cluster_col column in metadata storing cluster information
#' @param EnsDb.genome Ensembl annotation package of genome of interest (default = EnsDb.Hsapiens.v86)
#' @param blacklist_gr GenomicRanges object of ENCODE blacklisted regions for genome of interest
#' (default uses blacklist annotation implemented in Signac)
#'
#' @return Seurat object with qc metric stored in `meta.features` of `peaks` assay
#'
#' @import Seurat
#' @import Signac
#' @import dplyr
#' @import GenomicRanges
#' @importFrom tibble rownames_to_column
#' @importFrom tibble column_to_rownames
#' @importFrom matrixStats rowMaxs
#'
#' @export
compute_peakQC <- function(win_seu, cluster_col="seurat_clusters", EnsDb.genome = EnsDb.Hsapiens.v86, blacklist_gr = blacklist_hg38){
  ## Coverage stats
  peaks_count <- rowSums(win_seu@assays$peaks@counts)
  bin_peaks_mat <- win_seu@assays$peaks@counts
  bin_peaks_mat@x <- ifelse(bin_peaks_mat@x > 1, 1, bin_peaks_mat@x)
  peaks_feature_count <- rowSums(bin_peaks_mat)

  peaks_stats <- data.frame(tot_count=peaks_count, tot_cells = peaks_feature_count,
                            frac_cells=peaks_feature_count/ncol(bin_peaks_mat) ) %>%
    rownames_to_column("peak_id")

  ## Max fraction of accessible cells per cluster
  cl_frac <- fracAccessible_cluster(win_seu, assay="peaks", cluster_col)
  max_frac <- rowMaxs(cl_frac)
  names(max_frac) <- rownames(cl_frac)
  peaks_stats <- peaks_stats %>% mutate(max_frac = max_frac[peak_id])
  ## Annotations
  gr_peaks <- StringToGRanges(rownames(win_seu), sep=c(":","-"))
  gr_peaks <- annotate_gr(gr_peaks, EnsDb.genome, blacklist_gr)
  anno_df <- data.frame(gr_peaks@elementMetadata) %>%
    mutate(peak_id = rownames(win_seu))
  peaks_stats <- full_join(peaks_stats, anno_df, by="peak_id")
  win_seu@assays$peaks@meta.features <- peaks_stats
  win_seu@assays$peaks@meta.features <- win_seu@assays$peaks@meta.features %>% column_to_rownames("peak_id")
  return(win_seu)
}

#' Get annotation overlap
#'
#' Computes overlap between GenomicRanges object and annotation coordinates
#'
#' @param gr GenomicRanges object
#' @param anno_gr GenomicRanges object of annotations (usually loaded from ensembldb package)
#' @param anno_name character of name of annotation (e.g. genes, exons, promoters ...)
#' @param minoverlap minimum no of overlapping bps to call an overlap
#'
#' @return GenomicRanges object with overlap stored in metadata column
#'
get_annotation_overlap <- function(gr, anno_gr, anno_name, minoverlap = 10){
  ## Clean annotation object
  seqlevelsStyle(anno_gr) <- 'UCSC'
  anno_gr <- keepStandardChromosomes(anno_gr, pruning.mode = 'coarse')
  ## Find overlaps
  overlaps <- findOverlaps(gr, anno_gr, minoverlap = minoverlap)
  gr@elementMetadata[anno_name] <- 0
  gr@elementMetadata[anno_name][queryHits(overlaps),] <- 1
  return(gr)
}

#' Annotate GenomicRanges
#'
#' add width and overlap with annotations to metadata columns of GenomicRanges object
#'
#' @param gr GenomicRanges object to annotate
#' @param EnsDb.genome Ensembl annotation package to use for annotation. See details for installation (default = EnsDb.Hsapiens.v86, for hg38).
#' @param blacklist_gr GenomicRanges object of ENCODE blacklist regions (default: blacklist_hg38 from Signac package)
#' @param minoverlap minimum no of overlapping bps to call an overlap
#'
#' @details You can install EnsDb annotation packages by running
#' \code{BiocManager::install("EnsDb.Hsapiens.v86")}
#'
#' @import ensembldb
#' @import dplyr
#' @import GenomicRanges
#' @import GenomeInfoDb
#' @import Signac
#' @import ensembldb
#' @importFrom stats setNames
#'
#' @export
annotate_gr <- function(gr, EnsDb.genome = EnsDb.Hsapiens.v86, blacklist_gr = Signac::blacklist_hg38, minoverlap=10){
  ## Compute width
  gr$peak_width <- width(gr)

  ## Compute overlap with annotations
  exon_coords <- ensembldb::exons(EnsDb.genome, filter = ~ gene_biotype == "protein_coding")
  genes_coords <- ensembldb::genes(EnsDb.genome, filter = ~ gene_biotype == "protein_coding")
  promoter_coords <- ensembldb::promoters(EnsDb.genome, filter = ~ gene_biotype == "protein_coding", upstream=2000, downstream = 100)
  anno_list <- list(exon=exon_coords, gene=genes_coords, promoter=promoter_coords)
  for (i in seq_along(anno_list)) {
    gr <- get_annotation_overlap(gr, anno_gr = anno_list[[i]], anno_name = names(anno_list)[i], minoverlap=minoverlap)
  }
  gr$annotation <- data.frame(gr@elementMetadata) %>%
    mutate(annotation = case_when(exon==1 ~ "exon",
                                  promoter==1 ~ "promoter",
                                  gene==1 & exon==0 ~ "intron",
                                  TRUE ~ "intergenic"
                                  )) %>%
    pull(annotation)

  ## Overlapping gene
  seqlevelsStyle(genes_coords) <- 'UCSC'
  genes_coords <- keepStandardChromosomes(genes_coords, pruning.mode = 'coarse')
  genespromoters_coords <- Extend(genes_coords, upstream = 2000)
  overlap <- findOverlaps(gr, genespromoters_coords, minoverlap = minoverlap)
  gr$gene_name <- NA
  gr$gene_id <- NA
  gr$gene_name[queryHits(overlap)] <- genespromoters_coords$gene_name[subjectHits(overlap)]
  gr$gene_id[queryHits(overlap)] <- genespromoters_coords$gene_id[subjectHits(overlap)]

  ## Distance to nearest TSS
  tss_coords <- resize(unlist(range(split(genes_coords, ~ gene_id))), width=1)
  tss_distance <- distanceToNearest(gr, tss_coords)
  gr$tss_distance <- tss_distance@elementMetadata$distance

  ## Overlap with blacklist
  gr <- get_annotation_overlap(gr, anno_gr = blacklist_gr, anno_name = "ENCODE_blacklist")

  return(gr)
  }

#' Calculate accessible fraction per cluster
#'
#' @description For each feature it computes the fraction of cells per cluster in which that feature is accessible
#'
#' @param win_seu Seurat object
#' @param assay assay in seurat object to use
#' @param group_col meta.data column storing cluster identity (or other grouping variable for cells)
#'
#' @return matrix of features x clusters storing fractions
#'
#' @import ensembldb
#' @import dplyr
#' @import GenomicRanges
#' @import GenomeInfoDb
#' @import Signac
#' @import ensembldb
#'
#' @export
fracAccessible_cluster <- function(win_seu, assay, group_col="seurat_clusters"){
  ## Aggregate accessibility by cluster
  bin_seu <- BinarizeCounts(win_seu, assay)
  cl_sum <- sapply(unique(bin_seu@meta.data[,group_col]), function(x){
    cl_ix = which(bin_seu@meta.data[,group_col]==x)
    return(rowSums(GetAssay(bin_seu, assay)@counts[,cl_ix]))
  })
  colnames(cl_sum) <- paste0("cl_", unique(win_seu$seurat_clusters))
  ## Normalize by no of cells in cluster
  cl_tot <- table(bin_seu$seurat_clusters)
  cl_tot <- setNames(cl_tot, paste0("cl_",names(cl_tot)))
  cl_frac <- sapply(colnames(cl_sum), function(x) cl_sum[,x] / cl_tot[x])
  colnames(cl_frac) <- paste0("fracAccessible_", colnames(cl_frac))
  return(cl_frac)
}




