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
library(openxlsx)
library(writexl)
library(ggplotify)

source('/home/sjcho/yard/functions/R/save_ggplot2_to_ppt.R')
source('/home/sjcho/projects/utils.R')
source('/home/sjcho/projects/4-1BB/publication_figures/utils.R')
source('/home/sjcho/yard/functions/R/draw_contour.R')

save_path = '/home/sjcho/projects/4-1BB/publication_figures/new/figure3./outs_publish/CD8_CD62L'
so.cd8 <- readRDS('/home/sjcho/projects/4-1BB/publication_figures/new/figure1./outs_with_Sell/so.cd8.rds')
so.cd4 <- readRDS('/home/sjcho/projects/4-1BB/publication_figures/new/figure1./outs_with_Sell/so.cd4.rds')
so.full <- readRDS('/home/sjcho/projects/4-1BB/publication_figures/new/figure1./outs_with_Sell/so.full.rds')
so.full <- add_clone_info(so.full)

saveRDS(so.full, file = '/home/sjcho/projects/4-1BB/publication_figures/preprocessed_seurat_obj.rds')

############## 0. 4-1BB exp level and proportion

p <- VlnPlot(subset(so.full, hash.merged == 'B16'), features = 'Tnfrsf9', pt.size = 0, group.by = 'ann_level2', cols = ann_level2_cols, alpha = 0.5) & geom_boxplot(width = 0.2, fill = 'white')
p <- p + NoLegend() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
saveRDS(p, paste0(save_path, '/Tnfrsf9_expression_by_ann_level2.rds'))
ggsave(p, file = paste0(save_path, '/supple5_Tnfrsf9_expression_by_ann_level2.png'), width = 4, height = 4)

p <- ggplot(subset(so.cd8, hash.merged != 'MC38')@meta.data, aes(x = hash.merged, fill = ann_level2)) +
  geom_bar(position = "fill") +
  labs(y = "Proportion") +
  theme_classic() +
  scale_fill_manual(values = ann_level2_cols)
ggsave(p, file = paste0(save_path, '/ann_level2_proportion_barplot_CD8.png'), width = 5, height = 4)
saveRDS(p, paste0(save_path, '/ann_level2_proportion_barplot_CD8.rds'))

############## 1. CD8
############## 1.1. CD8 UMAP
ann_level2_cols_cd8 = ann_level2_cols[levels(so.cd8$ann_level2)]

### total UMAP
p <- DimPlot(so.cd8, group.by = 'ann_level2', cols = ann_level2_cols_cd8, pt.size = 0.5) + set_UMAP
ggsave(p, file = paste0(save_path, '/main3_umap_ann_level2_CD8.png'), width = 6.5, height = 5)

### clone size
so.cd8$Clone_size_group = as.vector(so.full@meta.data[Cells(so.cd8), 'Clone_size_group'])
so.cd8$Clone_size_group = factor(so.cd8$Clone_size_group, levels = so.full$Clone_size_group %>% levels())

p <- DimPlot(so.cd8, group.by = 'Clone_size_group', cols = colors_clone_size,  
            pt.size = 0.75, split.by = 'hash.merged', ncol = 2) & set_UMAP # , order = rev(levels(so.full$Clone_size_group))
p <- p &
  theme(
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    strip.text = element_blank()
  )
ggsave(p, file = paste0(save_path, '/main3_umap_clone_size_splitted_CD8.png'), width = 9, height = 10)

############## 1.2. CD8 expansion pattern

############# proportion barplot
so.cd8$clone_size_plot = ifelse(so.cd8$Clone_size_group == 'Not Detected', 'Not Detected', ifelse(so.cd8$Clone_size_group == 'Singleton', 'Singlet', 'Expanded'))
so.cd8$clone_size_plot = factor(so.cd8$clone_size_plot, levels = c('Not Detected', 'Singlet', 'Expanded'))
so.cd8$ann_level1 = as.vector(so.cd8$ann_level2)
so.cd8$ann_level1[so.cd8$ann_level2 %in% c('CD8.eff.', 'CD8.eff.prolif.')] = 'CD8.effector'
so.cd8$ann_level1[so.cd8$ann_level2 %in% c('CD8.exh.', 'CD8.exh.prolif.', 'Tpex.Sell+', 'Tpex.Sell-')] = 'CD8.exhausted'

celltype_ann_level1 = c('CD8.effector', 'CD8.exhausted')
so.cd8.plot <- subset(so.cd8, cells = Cells(so.cd8)[so.cd8$hash.merged != 'MC38' & so.cd8$ann_level1 %in% celltype_ann_level1])

manual_cols = c('CD8.effector' = '#d3a50d', 'CD8.exhausted' = '#a70909')
for (celltype in celltype_ann_level1) {
  so.cd8.tmp <- subset(so.cd8.plot, subset = (ann_level1 == celltype))
  current_color <- manual_cols[celltype]

  bar_data <- so.cd8.tmp@meta.data %>%
    group_by(hash.merged) %>%
    summarise(
      n_expanded = sum(clone_size_plot == "Expanded"),
      n_singlet = sum(clone_size_plot == "Singlet"),
      ratio = n_expanded / n_singlet
    )
  b16_val <- bar_data %>% filter(hash.merged == "B16") %>% pull(ratio)

  p <- ggplot(bar_data, aes(x = hash.merged, y = ratio, fill = hash.merged)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = hash.cols) +
    geom_hline(yintercept = b16_val, linetype = "dotted", color = "red") +
    scale_y_continuous(limits = c(0, NA)) +
    labs(x = NULL, y = NULL, title = celltype) +
    theme_classic() +
    NoLegend() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    ) +
    scale_x_discrete(labels = c("Control", "B7-H4x4-1BB", "ICB", "Combi"))

  p <- p + scale_x_discrete(labels = c("Control", "B7-H4x4-1BB", "ICB", "Combi"))
  p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

  ggsave(paste0(save_path, '/expanded_singlet_ratio_', celltype, '_barplot.png'), p, width = 4, height = 4)
  saveRDS(p, paste0(save_path, '/expanded_singlet_ratio_', celltype, '_barplot.rds'))
}

p.list = list()
for (celltype in c('CD8.eff.', 'CD8.eff.prolif.', 'CD8.exh.', 'CD8.exh.prolif.', 'Tpex.Sell+', 'Tpex.Sell-')) {
  so.cd8.tmp <- subset(so.cd8.plot, subset = (ann_level2 == celltype))
  current_color <- ann_level2_cols[celltype]

  bar_data <- so.cd8.tmp@meta.data %>%
    group_by(hash.merged) %>%
    summarise(
      n_expanded = sum(clone_size_plot == "Expanded"),
      n_singlet = sum(clone_size_plot == "Singlet"),
      ratio = n_expanded / n_singlet
    )
  b16_val <- bar_data %>% filter(hash.merged == "B16") %>% pull(ratio)

  p <- ggplot(bar_data, aes(x = hash.merged, y = ratio, fill = hash.merged)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = hash.cols) +
    geom_hline(yintercept = b16_val, linetype = "dotted", color = "red") +
    scale_y_continuous(limits = c(0, NA)) +
    labs(x = NULL, y = NULL, title = celltype) +
    theme_classic() +
    NoLegend() +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    ) +
    scale_x_discrete(labels = c("Control", "B7-H4x4-1BB", "ICB", "Combi"))

  p <- p + scale_x_discrete(labels = c("Control", "B7-H4x4-1BB", "ICB", "Combi"))
  p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  p.list[[celltype]] = p

  print(bar_data)
}
p2 <- wrap_plots(p.list, ncol = 6)
saveRDS(p.list, paste0(save_path, '/expanded_singlet_ratio_ann_level2.rds'))
ggsave(p2, file = paste0(save_path, '/expanded_singlet_ratio_ann_level2.png'), width = 12, height = 4)

############## 1.3. effector signature scores
so.cd8$Clone_size_group = as.vector(so.full@meta.data[Cells(so.cd8), 'Clone_size_group'])
so.cd8$Clone_size_group = factor(so.cd8$Clone_size_group, levels = so.full$Clone_size_group %>% levels())
so.cd8$clone_size_plot = ifelse(so.cd8$Clone_size_group == 'Not Detected', 'Not Detected', ifelse(so.cd8$Clone_size_group == 'Singleton', 'Singlet', 'Expanded'))
so.cd8$clone_size_plot = factor(so.cd8$clone_size_plot, levels = c('Not Detected', 'Singlet', 'Expanded'))
so.cd8$ann_level1 = as.vector(so.cd8$ann_level2)
so.cd8$ann_level1[so.cd8$ann_level2 %in% c('CD8.eff.', 'CD8.eff.prolif.')] = 'CD8.effector'
so.cd8$ann_level1[so.cd8$ann_level2 %in% c('CD8.exh.', 'CD8.exh.prolif.', 'Tpex.Sell+', 'Tpex.Sell-')] = 'CD8.exhausted'

celltype_ann_level1 = c('CD8.effector', 'CD8.exhausted')

source('/home/sjcho/yard/functions/R/draw_stacked_violin_box.R')

library(fgsea)
library(msigdbr)
library(GO.db)
library(org.Mm.eg.db)
h_gene_sets = msigdbr(species = "mouse", category = "H")

go_list = list(
  'GO:0044194' = 'cytolytic granule',
  'GO:0001906' = 'cell killing',
  'GO:0031343' = 'positive regulation of cell killing',
  'GO:0001909' = 'leukocyte mediated cytotoxicity',
  'GO:0001913' = 'T cell mediated cytotoxicity',
  'GO:0006952' = 'defense response'
)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

so.subset <- subset(so.full, cells = Cells(so.full)[(so.full$ann_level1 %in% celltype_ann_level1) & (so.full$hash.merged %in% c('B16', 'bsAb', 'ICB', 'Combi'))])
go_cols <- c()

go_genes = list()
for (go_term in names(go_list)) {
    genes <- AnnotationDbi::select(org.Mm.eg.db,
                                   keys = go_term,
                                   columns = c("GO", "SYMBOL"),
                                   keytype = "GOALL")
    go_genes[[go_list[[go_term]]]] <- unique(genes$SYMBOL)
}

for (go in names(go_genes)) {
  signature_genes <- go_genes[[go]]
  signature_genes <- signature_genes[signature_genes %in% rownames(so.subset)]
  so.subset <- AddModuleScore(so.subset, features = list(signature_genes), name = go)
  go_cols <- c(go_cols, paste0(go, "1"))
}

so.subset$hash_w_celltype <- paste0(so.subset$hash.merged, "_", so.subset$ann_level1)
levels <- c()
for (celltype in celltype_ann_level1) {
    for (condition in c('B16', 'bsAb', 'ICB', 'Combi')) {
        levels = c(levels, paste0(condition, "_", celltype))
    }
}
so.subset$hash_w_celltype <- factor(so.subset$hash_w_celltype, levels = levels)

score_data <- FetchData(so.subset, vars = c("hash_w_celltype", go_cols))
avg_mat <- score_data %>%
  group_by(hash_w_celltype) %>%
  summarise(across(everything(), mean)) %>%
  as.data.frame()
  
rownames(avg_mat) <- avg_mat$hash_w_celltype
avg_mat <- as.matrix(avg_mat[, -1])
colnames(avg_mat) <- names(go_genes)
plot_mat <- t(avg_mat)
plot_mat <- t(apply(plot_mat, 1, scale))
rownames(plot_mat) <- names(go_genes)
colnames(plot_mat) <- rownames(avg_mat)

plot_long <- as.data.frame(plot_mat) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "Group", values_to = "Z") %>%
  mutate(Group = factor(Group, levels = colnames(plot_mat)),
         Gene = factor(Gene, levels = rev(rownames(plot_mat))))
max_val <- max(abs(plot_long$Z))
p <- ggplot(plot_long, aes(Group, Gene, fill = Z)) +
  geom_tile(color = "black", linewidth = 0.2) +
  scale_fill_distiller(palette = "RdBu", direction = -1, limit = c(-max_val, max_val)) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(title = NULL, x = NULL, y = NULL, fill = "Z-score")
ggsave(p, file = paste0(save_path, '/main3_eff_exh_signature_score_heatmap.png'), width = 6, height = 4)
ggsave(p, file = paste0(save_path, '/main3_eff_exh_signature_score_heatmap.pdf'), width = 6, height = 4)
ggplot2pptx(p, file = paste0(save_path, '/main3_eff_exh_signature_score_heatmap.pptx'), width = 6, height = 4)

raw_data_for_excel <- p$data
write_xlsx(raw_data_for_excel, paste0(save_path, "/excel_main3_Heatmap_RawData_eff_exh.xlsx"))

############## 1.3.2 effector signature / ann_level_2
celltype_ann_level2 = c('CD8.eff.', 'CD8.eff.prolif.', 'CD8.exh.', 'CD8.exh.prolif.')

so.subset <- subset(so.full, cells = Cells(so.full)[(so.full$ann_level2 %in% celltype_ann_level2) & (so.full$hash.merged %in% c('B16', 'bsAb', 'ICB', 'Combi'))])
go_cols <- c()

for (go in names(go_genes)) {
  signature_genes <- go_genes[[go]]
  signature_genes <- signature_genes[signature_genes %in% rownames(so.subset)]
  so.subset <- AddModuleScore(so.subset, features = list(signature_genes), name = go)
  go_cols <- c(go_cols, paste0(go, "1"))
}

so.subset$hash_w_celltype <- paste0(so.subset$hash.merged, "_", so.subset$ann_level2)
levels <- c()
for (celltype in celltype_ann_level2) {
    for (condition in c('B16', 'bsAb', 'ICB', 'Combi')) {
        levels = c(levels, paste0(condition, "_", celltype))
    }
}
so.subset$hash_w_celltype <- factor(so.subset$hash_w_celltype, levels = levels)

score_data <- FetchData(so.subset, vars = c("hash_w_celltype", go_cols))
avg_mat <- score_data %>%
  group_by(hash_w_celltype) %>%
  summarise(across(everything(), mean)) %>%
  as.data.frame()
  
rownames(avg_mat) <- avg_mat$hash_w_celltype
avg_mat <- as.matrix(avg_mat[, -1])
colnames(avg_mat) <- names(go_genes)
plot_mat <- t(avg_mat)
plot_mat <- t(apply(plot_mat, 1, scale))
rownames(plot_mat) <- names(go_genes)
colnames(plot_mat) <- rownames(avg_mat)

plot_long <- as.data.frame(plot_mat) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(-Gene, names_to = "Group", values_to = "Z") %>%
  mutate(Group = factor(Group, levels = colnames(plot_mat)),
         Gene = factor(Gene, levels = rev(rownames(plot_mat))))
max_val <- max(abs(plot_long$Z))

cond_display <- rep(c('B16', 'B7-H4x4-1BB', 'ICB', 'Combi'), length(celltype_ann_level2))
names(cond_display) <- colnames(plot_mat)

n_cols <- ncol(plot_mat)
n_genes <- nrow(plot_mat)

vline_pos <- seq(4.5, n_cols - 0.5, by = 4)
ct_pos   <- seq(2.5, n_cols, by = 4)
ct_label <- celltype_ann_level2

p <- ggplot(plot_long, aes(Group, Gene, fill = Z)) +
  geom_tile(color = "black", linewidth = 0.2) +
  geom_vline(xintercept = vline_pos, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  annotate("text", x = ct_pos, y = n_genes + 0.7,
           label = ct_label, size = 3, hjust = 0.5, fontface = "bold") +
  scale_fill_distiller(palette = "RdBu", direction = -1, limit = c(-max_val, max_val)) +
  scale_x_discrete(labels = cond_display) +
  coord_cartesian(clip = "off") +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(angle = 90, hjust = 1),
    axis.ticks   = element_blank(),
    panel.grid   = element_blank(),
    plot.margin  = margin(t = 40, r = 10, b = 10, l = 10)
  ) +
  labs(x = NULL, y = NULL, fill = "Z-score")
ggsave(p, file = paste0(save_path, '/supple5_LongAnnotation_signature_score_heatmap.png'), width = 10, height = 4)
ggsave(p, file = paste0(save_path, '/supple5_longAnnotation_signature_score_heatmap.pdf'), width = 10, height = 4)
ggplot2pptx(p, file = paste0(save_path, '/supple5_longAnnotation_signature_score_heatmap.pptx'), width = 10, height = 4)
saveRDS(p, paste0(save_path, '/supple5_LongAnnotation_signature_score_heatmap.rds'))

raw_data_for_excel <- p$data
write_xlsx(raw_data_for_excel, paste0(save_path, "/excel_supple5_Heatmap_RawData_eff_exh.xlsx"))

############## 1.4. effector molecule expression
genes = c('Tbx21', 'Ifng', 'Tnf', 'Gzma', 'Gzmb', 'Prf1')
celltype_ann_level1 = c('CD8.effector', 'CD8.exhausted')

so.plot <- subset(so.full, hash.merged != 'MC38' & ann_level1 %in% celltype_ann_level1)
so.plot$celltype_condition = paste0(so.plot$ann_level1, '_', so.plot$hash.merged)

celltype_condition_levels = c()
for (celltype in celltype_ann_level1) {
    for (condition in c('B16', 'bsAb', 'ICB', 'Combi')) {
        celltype_condition_levels = c(celltype_condition_levels, paste0(celltype, '_', condition))
    }
}
so.plot$celltype_condition = factor(so.plot$celltype_condition, levels = celltype_condition_levels)

p2 <- DotPlot(so.plot, group.by = 'celltype_condition', features = genes) +
      theme_classic() +
      labs(x = 'celltype_condition', y = 'genes') +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
      scale_size(range = c(2, 7)) + 
      coord_flip()
p2 <- p2 + geom_hline(yintercept = c(4.5), linetype = "dashed", color = "black")
ggsave(paste0(save_path, '/effector_molecule_dotplot_by_condition_and_celltype.png'), p2, width = 7, height = 5)
ggsave(paste0(save_path, '/effector_molecule_dotplot_by_condition_and_celltype.pdf'), p2, width = 7, height = 5)
ggplot2pptx(p2, paste0(save_path, '/effector_molecule_dotplot_by_condition_and_celltype.pptx'), width = 7, height = 5)

plot_data <- p2$data
write.xlsx(plot_data, paste0(save_path, '/effector_molecule_dotplot_raw_data.xlsx'))

############## 1.4.2 effector molecule expression / ann_level2
genes = c('Tbx21', 'Ifng', 'Tnf', 'Gzma', 'Gzmb', 'Prf1')
celltype_ann_level2 = c('CD8.eff.', 'CD8.eff.prolif.', 'CD8.exh.', 'CD8.exh.prolif.')

so.plot <- subset(so.full, hash.merged != 'MC38' & ann_level2 %in% celltype_ann_level2)
so.plot$celltype_condition = paste0(so.plot$ann_level2, '_', so.plot$hash.merged)

celltype_condition_levels = c()
for (celltype in celltype_ann_level2) {
    for (condition in c('B16', 'bsAb', 'ICB', 'Combi')) {
        celltype_condition_levels = c(celltype_condition_levels, paste0(celltype, '_', condition))
    }
}
so.plot$celltype_condition = factor(so.plot$celltype_condition, levels = celltype_condition_levels)

n_genes <- length(genes)
ct_pos  <- seq(2.5, length(celltype_condition_levels) - 1.5, by = 4)  # 2.5, 6.5, 10.5 ...

p2 <- DotPlot(so.plot, group.by = 'celltype_condition', features = genes) +
      theme_classic() +
      labs(x = NULL, y = NULL) +
      theme(
        axis.text.x  = element_text(angle = 90, hjust = 1),
        plot.margin  = margin(t = 40, r = 10, b = 10, l = 10)
      ) +
      coord_flip(clip = "off") +
      scale_size(range = c(2, 7)) +
      geom_hline(yintercept = seq(4.5, 16.5, by = 4), linetype = "dashed", color = "black") +
      scale_y_discrete(labels = rep(c("Control", "B7-H4x4-1BB", "ICB", "Combi"), length(celltype_ann_level2))) +
      annotate("text",
               x     = n_genes + 0.8,
               y     = ct_pos,
               label = celltype_ann_level2,
               size  = 3, hjust = 0.5, fontface = "bold") + set_color_RdBu

ggsave(paste0(save_path, '/effector_molecule_dotplot_by_condition_and_celltype_LongAnnotation.png'), p2, width = 14, height = 5)
ggsave(paste0(save_path, '/effector_molecule_dotplot_by_condition_and_celltype_LongAnnotation.pdf'), p2, width = 14, height = 5)
ggplot2pptx(p2, paste0(save_path, '/effector_molecule_dotplot_by_condition_and_celltype_LongAnnotation.pptx'), width = 14, height = 5)
saveRDS(p2, paste0(save_path, '/effector_molecule_dotplot_by_condition_and_celltype_LongAnnotation.rds'))

plot_data <- p2$data
write.xlsx(plot_data, paste0(save_path, '/effector_molecule_dotplot_by_condition_and_celltype_LongAnnotation_raw_data.xlsx'))

############## 1.4.3 effector molecule expression in vlnplot / ann_level2
genes = c('Tbx21', 'Ifng', 'Tnf', 'Gzma', 'Gzmb', 'Prf1')
celltype_ann_level2 = c('CD8.eff.Prf1', 'CD8.eff.Gzmk', 'CD8.eff.prolif.', 'CD8.exh.', 'CD8.pex', 'CD8.exh.prolif.')

so.plot <- subset(so.full, hash.merged != 'MC38' & ann_level2 %in% celltype_ann_level2)
so.plot$celltype_condition = paste0(so.plot$ann_level2, '_', so.plot$hash.merged)

celltype_condition_levels = c()
cols = c()
for (celltype in celltype_ann_level2) {
    for (condition in c('B16', 'bsAb', 'ICB', 'Combi')) {
        celltype_condition_levels = c(celltype_condition_levels, paste0(celltype, '_', condition))
        cols = c(cols, ann_level2_cols[[celltype]])
    }
}
so.plot$celltype_condition = factor(so.plot$celltype_condition, levels = celltype_condition_levels)

source('/home/sjcho/yard/functions/R/draw_stacked_violin_box.R')
p <- StackedVlnPlot_box(so.plot, group.by = 'celltype_condition', features = c('Prf1', 'Ifng', 'Tnf'), cols = cols)
ggsave(paste0(save_path, '/StackedVlnPlot_Box_CD8_effector_genes.png'), plot = p, width = 5, height = 10)
ggsave(paste0(save_path, '/StackedVlnPlot_Box_CD8_effector_genesr.pdf'), plot = p, width = 5, height = 10)
ggplot2pptx(p, paste0(save_path, '/StackedVlnPlot_Box_CD8_effector_genes.pptx'), width = 5, height = 10)