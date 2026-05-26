# conda activate project_lung_exercise_R
### select only CD8 T cell and singlet (determined by Hashtag)
### pool1 : 3 sample, so use HTO demux
### pool2 : 2 sample, so use manual threshold (determined on log10(raw-hashtag count))

library(Seurat)
library(pheatmap)
library(magrittr)
library(ggplot2)
library(dplyr)
library(scRepertoire)
library(RColorBrewer)
library(egg)
library(harmony)
library(future)
library(Matrix)
library(scales)
library(tibble)
library(tidyr)
library(patchwork)


source('/home/sjcho/yard/functions/R/save_ggplot2_to_ppt.R')
source('/home/sjcho/projects/utils.R')
source('/home/sjcho/projects/4-1BB/publication_figures/utils.R')
source('/home/sjcho/yard/functions/R/draw_contour.R')

save_path = '/home/sjcho/projects/4-1BB/publication_figures/new/figure1./publish_outs'
so.cd8 <- readRDS('/home/sjcho/projects/4-1BB/publication_figures/new/figure1./outs_with_Sell/so.cd8.rds')
so.cd4 <- readRDS('/home/sjcho/projects/4-1BB/publication_figures/new/figure1./outs_with_Sell/so.cd4.rds')
so.full <- readRDS('/home/sjcho/projects/4-1BB/publication_figures/new/figure1./outs_with_Sell/so.full.rds')
so.full <- add_clone_info(so.full)

############## 1. Full UMAP
p <- DimPlot(so.full, group.by = 'ann_level2', cols = ann_level2_cols, pt.size = 0.5) + set_UMAP
ggsave(p, file = paste0(save_path, '/fig1_ann_level2_umap_full.png'), width = 7.5, height = 5)
                                  
so.cd4$ann_level2 = as.vector(so.full@meta.data[Cells(so.cd4), 'ann_level2'])
so.cd4$ann_level2 = factor(so.cd4$ann_level2, levels = c('CD4.naive.like', 'CD4.Th1.like', 
                                                        'fTreg', 'fTreg.prolif.',
                                                        'cTreg', 'cTreg.prolif.', 'Treg.LN.homing', 'Treg.IFN.response'))

so.cd8$ann_level2 = as.vector(so.full@meta.data[Cells(so.cd8), 'ann_level2'])
so.cd8$ann_level2 = factor(so.cd8$ann_level2, levels = c('CD8.naive.like', 
                                                        'CD8.eff.', 'CD8.eff.prolif.',
                                                        'CD8.exh.', 'CD8.exh.prolif.', 'Tpex.Sell+', 'Tpex.Sell-'))
################# 2. marker dotplot
################ 2.A. marker genes in CD4 T cell
marker.genes <- list(
  "CD4.naive"         = c('Lef1', 'Sell', 'Tcf7', 'Ccr7'),
  "CD4.Th1"           = c('Tbx21', 'Ifng', 'Tnf', 'Cxcr3', 'Gzmk'),
  "Treg"              = c('Foxp3'),
  "fTreg"      = c('Nkg7', 'Ccl3', 'Ccl4', 'Ccl5', 'Cst7'),
  "cTreg"       = c('Ikzf2', 'Il2ra', 'Nrp1', 'Dock2', 'Zeb1'),
  "Treg.LN.homing"    = c('Klf2', 'S1pr1', 'Areg', 'Ccr2'),
  "Treg.IFN.response" = c('Ifit1', 'Ifit2'),
  "proliferating"     = c('Top2a', 'Mki67')
)
p <- DotPlot(so.cd4, group.by = 'ann_level2', features = marker.genes, scale = TRUE,  dot.scale = 4, col.max = 1.5, col.min = -1.5)
p <- p + scale_color_gradientn(colors = rev(brewer.pal(n = 9, name = "RdBu"))) + RotatedAxis()
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p <- p + scale_size(range = c(1, 6))
ggsave(p, file = paste0(save_path, '/fig_marker_genes_dotplot_CD4.png'), width = 12, height = 6)
ggsave(p, file = paste0(save_path, '/fig_marker_genes_dotplot_CD4.pdf'), width = 12, height = 6)
ggplot2pptx(p, file = paste0(save_path, '/fig_marker_genes_dotplot_CD4.pptx'), width = 12, height = 6)
saveRDS(p, file = paste0(save_path, '/fig_marker_genes_dotplot_CD4.rds'))

################# 3.B. marker genes in CD8 T cell
marker.genes <- list(
  "CD8.naive"         = c('Lef1', 'Sell', 'Tcf7', 'Ccr7'),
  "CD8.eff."     = c('Gzmb', 'Prf1', 'Gzmk', 'Ifng', 'Nkg7', 'Cxcr6'),
  "CD8.exh."    = c('Tox', 'Havcr2', 'Entpd1', 'Pdcd1', 'Lag3', 'Ctla4'),
  "proliferating"     = c('Top2a', 'Mki67')
)
p <- DotPlot(so.cd8, group.by = 'ann_level2', features = marker.genes, scale = TRUE,  dot.scale = 4, col.max = 1.5, col.min = -1.5)
p <- p + scale_color_gradientn(colors = rev(brewer.pal(n = 9, name = "RdBu"))) + RotatedAxis()
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
p <- p + scale_size(range = c(1, 6))
ggsave(p, file = paste0(save_path, '/fig_marker_genes_dotplot_CD8.png'), width = 12, height = 6)
ggsave(p, file = paste0(save_path, '/fig_marker_genes_dotplot_CD8.pdf'), width = 12, height = 6)
ggplot2pptx(p, file = paste0(save_path, '/fig_marker_genes_dotplot_CD8.pptx'), width = 12, height = 6)
saveRDS(p, file = paste0(save_path, '/fig_marker_genes_dotplot_CD8.rds'))

source('/home/sjcho/yard/functions/R/draw_violin_split.R')
so.full$ann_level_for_B16_MC38_compare = as.vector(so.full$ann_level2)
so.full$ann_level_for_B16_MC38_compare = factor(so.full$ann_level2, levels = c(levels(so.cd8$ann_level2), levels(so.cd4$ann_level2)))
so.full$ann_level1_for_B16_MC38_compare = as.vector(so.full$ann_level1)
so.full$ann_level1_for_B16_MC38_compare = factor(so.full$ann_level1, levels = rev(c('Treg', 'CD4.Th1.like', 'CD4.naive.like', 'CD8.naive.like', 'CD8.exhausted', 'CD8.effector')))

so.control <- subset(so.full, cells = Cells(so.full)[so.full$hash.merged %in% c('MC38', 'B16')])
so.control$hash.merged <- factor(as.vector(so.control$hash.merged), levels = c('MC38', 'B16'))

p <- somi_violinplot_split(so.control, "Tnfrsf9", "ann_level_for_B16_MC38_compare", "hash.merged", col=c('#c4635d', '#73a5ba'), y.mul = 1, stat=FALSE)
ggsave(p, file = paste0(save_path, '/fig_Tnfrsf9_vln_by_celltype_long_annotation.png'), width = 8, height = 4)
ggsave(p, file = paste0(save_path, '/fig_Tnfrsf9_vln_by_celltype_long_annotation.pdf'), width = 8, height = 4)
ggplot2pptx(p, file = paste0(save_path, '/fig_Tnfrsf9_vln_by_celltype_long_annotation.pptx'), width = 8, height = 4)

p <- somi_violinplot_split(so.control, "Tnfrsf9", "ann_level1_for_B16_MC38_compare", "hash.merged", col=c('#c4635d', '#73a5ba'), y.mul = 1, stat=FALSE)
ggsave(p, file = paste0(save_path, '/fig_Tnfrsf9_vln_by_celltype_ann1.png'), width = 6, height = 4)
ggsave(p, file = paste0(save_path, '/fig_Tnfrsf9_vln_by_celltype_ann1.pdf'), width = 6, height = 4)
ggplot2pptx(p, file = paste0(save_path, '/fig_Tnfrsf9_vln_by_celltype_ann1.pptx'), width = 6, height = 4)

saveRDS(so.full, file = paste0(save_path, '/so.full.rds'))
saveRDS(so.cd4, file = paste0(save_path, '/so.cd4.rds'))
saveRDS(so.cd8, file = paste0(save_path, '/so.cd8.rds'))