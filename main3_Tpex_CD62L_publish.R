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
library(ggsignif)
library(openxlsx)
library(reshape2)
library(ggplotify)

source('/home/sjcho/yard/functions/R/save_ggplot2_to_ppt.R')
source('/home/sjcho/projects/utils.R')
source('/home/sjcho/projects/4-1BB/publication_figures/utils.R')
source('/home/sjcho/yard/functions/R/draw_contour.R')
source('/home/sjcho/yard/functions/R/draw_volcano.R')

save_path = '/home/sjcho/projects/4-1BB/publication_figures/new/figure3./outs_publish/Tpex_CD62L'
so.cd8 <- readRDS('/home/sjcho/projects/4-1BB/publication_figures/new/figure1./outs_with_Sell/so.cd8.rds')
so.cd4 <- readRDS('/home/sjcho/projects/4-1BB/publication_figures/new/figure1./outs_with_Sell/so.cd4.rds')
so.full <- readRDS('/home/sjcho/projects/4-1BB/publication_figures/new/figure1./outs_with_Sell/so.full.rds')
so.full <- add_clone_info(so.full)

############## 1. subset, normalize, harmony ##############
so.Tpex = subset(so.cd8, ann_level2 %in% c('Tpex.Sell+', 'Tpex.Sell-'))
so.Tpex <- log_normalize(so.Tpex, save_path = save_path, data_name = 'Tpex', pcs = 7, nfeatures = 2000)
so.Tpex <- harmony_on_RNA(so.Tpex, batch = 'pool', pc.rna = 7, save_path = save_path, sample_name = 'Tpex')
so.Tpex <- FindClusters(so.Tpex, resolution = 2.5)

p <- DimPlot(so.Tpex, label = T, reduction = 'umap.harmony.rna') + set_UMAP  
ggsave(p, file = paste0(save_path, '/umap_Tpex_clusters_before_annotation.png'), width = 5, height = 4)

p <- FeaturePlot(so.Tpex, features = c('Sell'), reduction = 'umap.harmony.rna', ncol = 1, order = T, pt.size = 1) & set_UMAP_featureplot
p <- p + ggtitle('')
ggsave(p + NoLegend(), file = paste0(save_path, '/umap_Tpex_Sell.png'), width = 4, height = 4)

############### 2. visualize
############### 2.1. UMAP (condition, contour)
cols = c('#9d0200', '#7bffa0', '#7c36e4', '#a39600', '#004cba', '#e87700', '#0086b8', '#ffd463', '#9e0096', '#02a290', '#e956bd', '#fbe5ff', '#002612', '#ffbe8f', '#392200')
p <- DimPlot(so.Tpex, label = FALSE, reduction = 'umap.harmony.rna', cols = cols) + set_UMAP  
ggsave(p, file = paste0(save_path, '/umap_Tpex_clusters_before_annotation.png'), width = 5, height = 4)
saveRDS(p, file = paste0(save_path, '/umap_Tpex_clusters_before_annotation.rds'))
p <- DimPlot(so.Tpex, label = TRUE, reduction = 'umap.harmony.rna', cols = cols) + set_UMAP  
ggsave(p, file = paste0(save_path, '/umap_Tpex_clusters_before_annotation_label.png'), width = 5, height = 4)

p <- VlnPlot(so.Tpex, features = c('Sell', 'Myb'), group.by = 'seurat_clusters', pt.size = 0, cols = cols) & geom_boxplot(width = 0.2, alpha = 0.5)
ggsave(p, file = paste0(save_path, '/vln_Tpex_Myb_Sell_by_cluster.png'), width = 8, height = 4)

p <- fp_sjcho(so.Tpex, features = c('Sell'), reduction = 'umap.harmony.rna', ncol = 1, order = T, pt.size = 1) & set_UMAP_featureplot
ggsave(p, file = paste0(save_path, '/umap_Tpex_Sell.png'), width = 5, height = 4)
legend <- get_legend(p)
saveRDS(as.ggplot(legend), file = paste0(save_path, '/umap_Tpex_Sell_legend.rds'))
saveRDS(p, file = paste0(save_path, '/umap_Tpex_Sell.rds'))

so.Tpex@meta.data[, 'pool'] <- factor(so.Tpex@meta.data[ ,'pool'], levels = c('pool1', 'pool2'))
tmp_cols = rep('black', 2)
names(tmp_cols) = c('pool1', 'pool2')
p <- contour_plot_merge(so.Tpex, name_of_dimension = 'umap.harmony.rna', condition_name = 'pool', 
                  celltype_name = 'hash.merged', celltype_color = hash.cols,
                  condition_color = tmp_cols, nbin = 5, linewidth = 0.5, line_alpha = 0.75, alpha = 1, ncol = 2, nrow = 1)
ggsave(p, file = paste0(save_path, '/umap_contour_hash.merged_pool.png'), width = 5, height = 2.5)
ggsave(p, file = paste0(save_path, '/umap_contour_hash.merged_pool.pdf'), width = 5, height = 2.5)
ggplot2pptx(p, width = 8, height = 4, paste0(save_path, '/umap_contour_hash.merged_pool.pptx'))

tmp_cols = rep('black', 5)
names(tmp_cols) = levels(so.Tpex$'hash.merged')
p <- contour_plot_merge(so.Tpex, name_of_dimension = 'umap.harmony.rna', condition_name = 'hash.merged', 
                  celltype_name = 'hash.merged', celltype_color = hash.cols,
                  condition_color = tmp_cols, nbin = 5, linewidth = 0.5, line_alpha = 0.75, alpha = 1, ncol = 3, nrow = 2)
ggsave(p, file = paste0(save_path, '/umap_contour_hash.merged.png'), width = 7.5, height = 5)
ggsave(p, file = paste0(save_path, '/umap_contour_hash.merged.pdf'), width = 7.5, height = 5)
ggplot2pptx(p, width = 9, height = 6, paste0(save_path, '/umap_contour_hash.merged.pptx'))

so.Tpex$Tpex_subset = ifelse(so.Tpex$seurat_clusters %in% c(1, 9, 10, 3, 8, 13), 'Sell+ Tpex', 'Sell- Tpex')
so.Tpex$Tpex_subset = factor(so.Tpex$Tpex_subset, levels = c('Sell- Tpex', 'Sell+ Tpex'))

#################### 3. proportion
metadata <- so.Tpex@meta.data
p <- ggplot(metadata, aes(x = hash.merged, fill = Tpex_subset)) +
  geom_bar(position = "fill") +
  theme_classic() + 
  scale_fill_manual(values = c("Sell- Tpex" = "#c99090", "Sell+ Tpex" = "#370202")) + 
  scale_y_continuous(labels = scales::percent) + 
  labs(
    title = "Proportion of Tpex_subset by hash.merged",
    x = "Group (hash.merged)",
    y = "Proportion (%)",
    fill = "Tpex Subset"
  ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(p, file = paste0(save_path, '/proportion_Tpex_subset_by_hash.merged.png'), width = 4, height = 3)
ggsave(p, file = paste0(save_path, '/proportion_Tpex_subset_by_hash.merged.pdf'), width = 4, height = 3)
saveRDS(p, file = paste0(save_path, '/proportion_Tpex_subset_by_hash.merged.rds'))
ggplot2pptx(p, width = 5, height = 4, paste0(save_path, '/proportion_Tpex_subset_by_hash.merged.pptx'))

summary_data <- so.Tpex@meta.data %>%
  group_by(hash.merged, Tpex_subset) %>%
  summarise(n = n(), .groups = 'drop') %>%
  group_by(hash.merged) %>%
  mutate(proportion = n / sum(n))
write.xlsx(summary_data, file = paste0(save_path, '/excel_proportion_Tpex_subset_by_hash.merged.xlsx'))

################ 4. Stemness and exhaustion score
so.Tpex <- subset(so.cd8, ann_level2 %in% c('Tpex.Sell+', 'Tpex.Sell-'))
so.Tpex$ann_level2 <- as.vector(so.Tpex$ann_level2)
so.Tpex$ann_level2 = factor(so.Tpex$ann_level2, levels = c('Tpex.Sell+', 'Tpex.Sell-'))

p1 <- VlnPlot(so.Tpex, group.by = 'ann_level2', features = 'Tcf7', pt.size = 0, cols = c('#370202', '#c99090')) + geom_boxplot(width = 0.2, fill = 'white', alpha = 0.75, outlier.size = 0.5)
p1 <- p1 + ylim(0, 1.1 * max(so.Tpex@assays$RNA$data['Tcf7', ])) + theme(axis.title.y = element_blank()) + NoLegend()

p2 <- VlnPlot(so.Tpex, group.by = 'ann_level2', features = 'Sell', pt.size = 0, cols = c('#370202', '#c99090')) + geom_boxplot(width = 0.2, fill = 'white', alpha = 0.75, outlier.size = 0.5)
p2 <- p2 + ylim(0, 1.1 * max(so.Tpex@assays$RNA$data['Sell', ])) + theme(axis.title.y = element_blank()) + NoLegend()

p3 <- VlnPlot(so.Tpex, group.by = 'ann_level2', features = 'Myb', pt.size = 0, cols = c('#370202', '#c99090')) + geom_boxplot(width = 0.2, fill = 'white', alpha = 0.75, outlier.size = 0.5)
p3 <- p3 + ylim(0, 1.1 * max(so.Tpex@assays$RNA$data['Myb', ])) + theme(axis.title.y = element_blank()) + NoLegend()

p4 <- VlnPlot(so.Tpex, group.by = 'ann_level2', features = 'Tox', pt.size = 0, cols = c('#370202', '#c99090')) + geom_boxplot(width = 0.2, fill = 'white', alpha = 0.75, outlier.size = 0.5)
p4 <- p4 + ylim(0, 1.1 * max(so.Tpex@assays$RNA$data['Tox', ])) + theme(axis.title.y = element_blank()) + NoLegend()

p5 <- VlnPlot(so.Tpex, group.by = 'ann_level2', features = 'Havcr2', pt.size = 0, cols = c('#370202', '#c99090')) + geom_boxplot(width = 0.2, fill = 'white', alpha = 0.75, outlier.size = 0.5)
p5 <- p5 + ylim(0, 1.1 * max(so.Tpex@assays$RNA$data['Havcr2', ])) + theme(axis.title.y = element_blank()) + NoLegend()

p6 <- VlnPlot(so.Tpex, group.by = 'ann_level2', features = 'Ctla4', pt.size = 0, cols = c('#370202', '#c99090')) + geom_boxplot(width = 0.2, fill = 'white', alpha = 0.75, outlier.size = 0.5)
p6 <- p6 + ylim(0, 1.1 * max(so.Tpex@assays$RNA$data['Ctla4', ])) + theme(axis.title.y = element_blank()) + NoLegend()
p <- cowplot::plot_grid(p1, p2, p3, p4, p5, p6, nrow = 1, ncol = 6)
saveRDS(p, file = paste0(save_path, '/vln_Tpex_stemness_exh_by_Tpex_subset.rds'))
saveRDS(p1, file = paste0(save_path, '/00vln_Tpex_Tcf7_by_Tpex_subset.rds'))
saveRDS(p2, file = paste0(save_path, '/00vln_Tpex_Sell_by_Tpex_subset.rds'))
saveRDS(p3, file = paste0(save_path, '/00vln_Tpex_Myb_by_Tpex_subset.rds'))
saveRDS(p4, file = paste0(save_path, '/00vln_Tpex_Tox_by_Tpex_subset.rds'))
saveRDS(p5, file = paste0(save_path, '/00vln_Tpex_Havcr2_by_Tpex_subset.rds'))
saveRDS(p6, file = paste0(save_path, '/00vln_Tpex_Ctla4_by_Tpex_subset.rds'))

p <- DotPlot(so.Tpex, features = c('Tcf7', 'Sell', 'Il7r', 'Myb', 'Tox', 'Havcr2', 'Ctla4'), group.by = 'Tpex_subset') + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p <- p + set_color_RdBu + scale_size(range = c(2, 6))
ggsave(p, file = paste0(save_path, '/dotplot_Tpex_stemness_exh_by_Tpex_subset.png'), width = 6, height = 4)

##################### 4.2. Sell+ Tpex
p <- DotPlot(subset(so.Tpex, Tpex_subset == 'Sell+ Tpex'), features = c('Tcf7', 'Sell', 'Il7r', 'Myb', 'Tox', 'Havcr2', 'Ctla4'), group.by = 'hash.merged') + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p <- p + set_color_RdBu + scale_size(range = c(2, 6))
ggsave(p, file = paste0(save_path, '/dotplot_Tpex_stemness_exh_Sell_p_Tpex_by_cluster.png'), width = 6, height = 4)

##################### stemness in naive-like
p <- DotPlot(subset(so.cd8, ann_level2 == 'CD8.naive.like' & hash.merged != 'MC38'), features = c('Tcf7', 'Myb', 'Lef1', 'Sell', 'Bcl6'), group.by = 'hash.merged')
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + set_color_RdBu
p <- p + coord_flip()
saveRDS(p, file = paste0(save_path, '/dotplot_cd8_naivelike_stemness.rds'))