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

source('/home/sjcho/yard/functions/R/save_ggplot2_to_ppt.R')
source('/home/sjcho/projects/utils.R')
source('/home/sjcho/projects/4-1BB/publication_figures/utils.R')
source('/home/sjcho/yard/functions/R/draw_contour.R')

save_path = '/home/sjcho/projects/4-1BB/publication_figures/new/figure3./outs_publish/CD4'
so.cd8 <- readRDS('/home/sjcho/projects/4-1BB/publication_figures/new/figure1./outs_with_Sell/so.cd8.rds')
so.cd4 <- readRDS('/home/sjcho/projects/4-1BB/publication_figures/new/figure1./outs_with_Sell/so.cd4.rds')
so.full <- readRDS('/home/sjcho/projects/4-1BB/publication_figures/new/figure1./outs_with_Sell/so.full.rds')
so.full <- add_clone_info(so.full)

############## 1. CD4
############## 1.1. CD4 UMAP
### proportion
p <- ggplot(subset(so.cd4, hash.merged != 'MC38')@meta.data, aes(x = hash.merged, fill = ann_level2)) +
  geom_bar(position = "fill") +
  labs(y = "Proportion") +
  theme_classic() +
  scale_fill_manual(values = ann_level2_cols)
ggsave(p, file = paste0(save_path, '/ann_level2_proportion_barplot_CD4.png'), width = 5, height = 4)
saveRDS(p, paste0(save_path, '/ann_level2_proportion_barplot_CD4.rds'))

### total UMAP
p <- DimPlot(so.cd4, group.by = 'ann_level2', cols = ann_level2_cols_cd4, pt.size = 0.5) + set_UMAP
ggsave(p, file = paste0(save_path, '/ann_level2_umap_CD4.png'), width = 7, height = 5)

### clone size
so.cd4$Clone_size_group = as.vector(so.full@meta.data[Cells(so.cd4), 'Clone_size_group'])
so.cd4$Clone_size_group = factor(so.cd4$Clone_size_group, levels = so.full$Clone_size_group %>% levels())

p <- DimPlot(so.cd4, group.by = 'Clone_size_group', cols = colors_clone_size, 
            pt.size = 0.75, split.by = 'hash.merged', ncol = 3, order = rev(levels(so.full$Clone_size_group))) & set_UMAP # 
p <- p &
  labs(title = NULL, subtitle = NULL) &
  theme(
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    strip.text = element_blank()
  )
ggsave(p, file = paste0(save_path, '/clone_size_splitted_umap_CD4.png'), width = 10, height = 5)

############# proportion barplot
so.cd4$Clone_size_group = as.vector(so.full@meta.data[Cells(so.cd4), 'Clone_size_group'])
so.cd4$Clone_size_group = factor(so.cd4$Clone_size_group, levels = so.full$Clone_size_group %>% levels())

so.cd4$clone_size_plot = ifelse(so.cd4$Clone_size_group == 'Not Detected', 'Not Detected', ifelse(so.cd4$Clone_size_group == 'Singleton', 'Singlet', 'Expanded'))
so.cd4$clone_size_plot = factor(so.cd4$clone_size_plot, levels = c('Not Detected', 'Singlet', 'Expanded'))
so.cd4$ann_level1 = as.vector(so.cd4$ann_level2)
so.cd4$ann_level1[so.cd4$ann_level2 %in% c('cTreg', 'cTreg.prolif.')] = 'Treg.stable'
so.cd4$ann_level1[so.cd4$ann_level2 %in% c('fTreg', 'fTreg.prolif.')] = 'Treg.fragile'

celltypes = c('Treg.stable', 'Treg.fragile', 'CD4.Th1.like')
so.cd4.plot <- subset(so.cd4, cells = Cells(so.cd4)[so.cd4$hash.merged != 'MC38' & so.cd4$ann_level1 %in% celltypes])
manual_cols = c('Treg.stable' = '#15574d', 'Treg.fragile' = '#00a0fc')

so.cd4.ftreg <- subset(so.cd4.plot, subset = (ann_level1 == "Treg.fragile"))
ftreg_max <- so.cd4.ftreg@meta.data %>%
  group_by(hash.merged) %>%
  summarise(ratio = sum(clone_size_plot == "Expanded") / sum(clone_size_plot == "Singlet")) %>%
  pull(ratio) %>% max()

for (celltype in celltypes) {
  so.cd4.tmp <- subset(so.cd4.plot, subset = (ann_level1 == celltype))
  current_color <- manual_cols[celltype]

  bar_data <- so.cd4.tmp@meta.data %>%
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
    scale_y_continuous(limits = c(0, ftreg_max * 1.05)) +
    labs(
      title = paste0("Expanded/Singlet Ratio: ", celltype),
      x = NULL, y = "Expanded / Singlet ratio"
    ) +
    theme_classic() +
    NoLegend() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    scale_x_discrete(labels = c("Control", "B7-H4x4-1BB", "ICB", "Combi"))

  ggsave(paste0(save_path, '/expanded_singlet_ratio_', celltype, '_barplot.png'), p, width = 4, height = 4)
  saveRDS(p, paste0(save_path, '/expanded_singlet_ratio_', celltype, '_barplot.rds'))
}

############## total proportion
df <- subset(so.cd4, cells = Cells(so.cd4)[!(so.cd4$hash.merged == 'MC38')])@meta.data %>%
  count(hash.merged, ann_level2) %>%
  group_by(hash.merged) %>%
  mutate(prop = n / sum(n))

p <- ggplot(df, aes(x = hash.merged, y = prop, fill = hash.merged)) +
  geom_col() +
  scale_fill_manual(values = hash.cols) +
  facet_wrap(~ann_level2, ncol = 4, scales = "free") +
  labs(y = "Proportion") +
  theme_classic() +
  guides(fill = "none") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
saveRDS(p, paste0(save_path, '/proportion_ann_level2_by_condition.rds'))

### assemble by patchwork
p1 <- wrap_plots(p.list, ncol = 7)

ggsave(paste0(save_path, '/relative_expansion_ratio_by_condition_expansion_over2.ann2.png'), p1, width = 21, height = 3)
ggsave(paste0(save_path, '/relative_expansion_ratio_by_condition_expansion_over2.ann2.pdf'), p1, width = 21, height = 3)
ggplot2pptx(p1, paste0(save_path, '/relative_expansion_ratio_by_condition_expansion_over2.ann2.pptx'), width = 21, height = 3)
write.csv(result_df, file = paste0(save_path, '/relative_expansion_ratio_by_condition_expansion_over2.ann2.csv'))

so.cd4.plot <- subset(so.cd4, cells = Cells(so.cd4)[so.cd4$hash.merged != 'MC38' & so.cd4$ann_level2 %in% c('fTreg', 'fTreg.prolif.', 'cTreg', 'cTreg.prolif.')])
so.cd4.plot$ann_level1 = as.vector(so.cd4.plot$ann_level2)
so.cd4.plot$ann_level1[so.cd4.plot$ann_level1 %in% c('fTreg', 'fTreg.prolif.')] = 'fTreg'
so.cd4.plot$ann_level1[so.cd4.plot$ann_level1 %in% c('cTreg', 'cTreg.prolif.')] = 'cTreg'

p.list = list()
for (celltype in c('fTreg', 'cTreg')) {
  so.cd4.tmp <- subset(so.cd4.plot, subset = (ann_level1 == celltype & hash.merged != 'MC38'))
  bar_data <- so.cd4.tmp@meta.data %>%
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
p2 <- wrap_plots(p.list, ncol = 2)
saveRDS(p.list, paste0(save_path, '/expanded_singlet_ratio_ann_level1_cd4.rds'))
ggsave(paste0(save_path, '/expanded_singlet_ratio_ann_level1_cd4.png'), p2, width = 21, height = 3)

############## 1.3. effector signature score
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

go_genes = list()
for (go_term in names(go_list)) {
    genes <- AnnotationDbi::select(org.Mm.eg.db,
                                   keys = go_term,
                                   columns = c("GO", "SYMBOL"),
                                   keytype = "GOALL")
    go_genes[[go_list[[go_term]]]] <- unique(genes$SYMBOL)
}

celltype_ann_level2 = c('fTreg', 'fTreg.prolif.')
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
ggsave(p, file = paste0(save_path, '/main4_fTreg_signature_score_heatmap.png'), width = 6, height = 4)
ggsave(p, file = paste0(save_path, '/main4_fTreg_signature_score_heatmap.pdf'), width = 6, height = 4)
ggplot2pptx(p, file = paste0(save_path, '/main4_fTreg_signature_score_heatmap.pptx'), width = 6, height = 4)
raw_data_for_excel <- p$data
saveRDS(p, paste0(save_path, "/main4_fTreg_signature_score_heatmap.rds"))
write_xlsx(raw_data_for_excel, paste0(save_path, "/excel_main4_fTreg_Heatmap_RawData_eff_exh.xlsx"))

############## 1.3.2 effector signature / CD4.Th1.like
celltype_ann_level2 = c('CD4.Th1.like')
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
ggsave(p, file = paste0(save_path, '/supple7_signature_score.Th1.like.png'), width = 10, height = 4)
ggsave(p, file = paste0(save_path, '/supple7_signature_score.Th1.like.pdf'), width = 10, height = 4)
saveRDS(p, paste0(save_path, "/supple7_signature_score.Th1.like.rds"))

raw_data_for_excel <- p$data
write_xlsx(raw_data_for_excel, paste0(save_path, "/excel_supple7_Heatmap_RawData_eff_exh.xlsx"))

############## 1.4. effector molecule expression
genes = c('Tbx21', 'Ifng', 'Tnf', 'Nkg7', 'Gzmb', 'Prf1')
# celltype_ann_level1 = c('fTreg', 'fTreg.prolif.')
celltype_ann_level1 = c('CD4.Th1.like')

so.plot <- subset(so.cd4, hash.merged != 'MC38' & ann_level2 %in% celltype_ann_level1)
p2 <- DotPlot(so.plot, group.by = 'hash.merged', features = genes) +
      theme_classic() +
      labs(x = 'celltype_condition', y = 'genes') +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
      scale_size(range = c(1, 6)) +
      coord_flip()
ggsave(paste0(save_path, '/main4_CD4.Th1_effector_molecule_dotplot_by_condition.png'), p2, width = 4, height = 5)
ggsave(paste0(save_path, '/main4_CD4.Th1_effector_molecule_dotplot_by_condition.pdf'), p2, width = 4, height = 5)
saveRDS(p2, paste0(save_path, '/main4_CD4.Th1_effector_molecule_dotplot_by_condition.rds'))
ggplot2pptx(p2, paste0(save_path, '/main4_fTreg_effector_molecule_dotplot_by_condition.pptx'), width = 4, height = 5)

source('/home/sjcho/yard/functions/R/draw_stacked_violin_box.R')
p <- StackedVlnPlot_box(so.plot, group.by = 'hash.merged', features = genes, cols = hash.cols)
ggsave(paste0(save_path, '/main4_fTreg_effector_molecule_StackedVlnPlot_Box_by_condition.png'), plot = p, width = 5, height = 10)
ggsave(paste0(save_path, '/main4_fTreg_effector_molecule_StackedVlnPlot_Box_by_condition.pdf'), plot = p, width = 5, height = 10)
ggplot2pptx(p, paste0(save_path, '/main4_fTreg_effector_molecule_StackedVlnPlot_Box_by_condition.pptx'), width = 5, height = 10)

so.plot <- subset(so.cd4, hash.merged != 'MC38' & ann_level2 %in% c('fTreg', 'CD4.Th1.like'))
so.plot$celltype_condition = paste0(so.plot$ann_level2, '_', so.plot$hash.merged)
celltype_condition_levels = c()
for (celltype in c('fTreg', 'CD4.Th1.like')) {
    for (condition in c('B16', 'bsAb', 'ICB', 'Combi')) {
        celltype_condition_levels = c(celltype_condition_levels, paste0(celltype, '_', condition))
    }
}
so.plot$celltype_condition = factor(so.plot$celltype_condition, levels = celltype_condition_levels)
p <- DotPlot(so.plot, group.by = 'celltype_condition', features = genes)
p <- p + theme_classic() +
      labs(x = 'celltype_condition', y = 'genes') +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
      coord_flip() +
      scale_size(range = c(1, 6))


############## 1.4.3 effector molecule expression / CD4.Th1.like
genes = c('Tbx21', 'Ifng', 'Tnf', 'Nkg7', 'Gzmb', 'Prf1')
so.plot <- subset(so.cd4, (hash.merged != 'MC38') & (ann_level2 %in% c('CD4.Th1.like')))
so.plot$celltype_condition = paste0(so.plot$ann_level2, '_', so.plot$hash.merged)

celltype_condition_levels = c()
for (celltype in celltype_ann_level2) {
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
      coord_flip() +
      scale_y_discrete(labels = rep(c("B16", "bsAb", "ICB", "Combi"), length(celltype_ann_level2))) +
      scale_size(range = c(1, 6))

ggsave(paste0(save_path, '/effector_molecule_dotplot_by_condition_Th1_like.png'), p2, width = 5, height = 5)
ggsave(paste0(save_path, '/effector_molecule_dotplot_by_condition_Th1_like.pdf'), p2, width = 5, height = 5)
saveRDS(p2, paste0(save_path, '/effector_molecule_dotplot_by_condition_Th1_like.rds'))

ggplot2pptx(p2, paste0(save_path, '/effector_molecule_dotplot_by_condition_and_celltype_ann_level2.pptx'), width = 5, height = 5)

data <- p2$data
write_xlsx(data, paste0(save_path, "/excel_effector_dotplot_RawData_ann_level2.xlsx"))

############## 1.4.4 effector molecule expression in vlnplot / ann_level2
celltype_ann_level2 = c('CD4.Th1.like', 'Treg.fragile', 'Treg.fragile.prolif.')

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
ggsave(paste0(save_path, '/StackedVlnPlot_Box_CD4_effector_genes.png'), plot = p, width = 5, height = 10)
ggsave(paste0(save_path, '/StackedVlnPlot_Box_CD4_effector_genesr.pdf'), plot = p, width = 5, height = 10)
ggplot2pptx(p, paste0(save_path, '/StackedVlnPlot_Box_CD4_effector_genes.pptx'), width = 5, height = 10)

############## 1.5. suppressive molecule expression
genes = c( 'Ikzf2', 'Il2ra', 'Ctla4', 'Tnfrsf18', 'Il10')
celltype_ann_level2 = c('cTreg', 'cTreg.prolif.', 'fTreg', 'fTreg.prolif.', 'Treg.LN.homing', 'Treg.IFN.response')

so.plot <- subset(so.full, hash.merged != 'MC38' & ann_level2 %in% c('cTreg', 'cTreg.prolif.'))
p2 <- DotPlot(so.plot, group.by = 'hash.merged', features = genes) +
      theme_classic() +
      labs(x = 'hash.merged', y = 'genes') +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
      coord_flip() +
      scale_size(range = c(1.5, 6))
p2 <- p2 + scale_y_discrete(labels = c("B16", "B7-H4x4-1BB", "ICB", "Combi"))
ggsave(paste0(save_path, '/suppressive_dotplot_by_condition_and_celltype.png'), p2, width = 4, height = 5)
ggsave(paste0(save_path, '/suppressive_dotplot_by_condition_and_celltype.pdf'), p2, width = 4, height = 5)
saveRDS(p2, paste0(save_path, '/suppressive_dotplot_by_condition_and_celltype.rds'))
ggplot2pptx(p2, paste0(save_path, '/suppressive_dotplot_by_condition_and_celltype.pptx'), width = 4, height = 5)

data <- p2$data
write_xlsx(data, paste0(save_path, "/excel_suppressive_dotplot_RawData.xlsx"))

############## 1.5.2 suppressive molecule expression / ann_level2
genes = c('Ikzf2', 'Il2ra', 'Ctla4', 'Tnfrsf18', 'Il10')
celltype_ann_level2 = c('cTreg', 'cTreg.prolif.', 'fTreg', 'fTreg.prolif.', 'Treg.LN.homing', 'Treg.IFN.response')

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
ct_pos  <- seq(2.5, length(celltype_condition_levels) - 1.5, by = 4)

p2 <- DotPlot(so.plot, group.by = 'celltype_condition', features = genes) +
      theme_classic() +
      labs(x = NULL, y = NULL) +
      theme(
        axis.text.x  = element_text(angle = 90, hjust = 1),
        plot.margin  = margin(t = 40, r = 10, b = 10, l = 10)
      ) +
      coord_flip(clip = "off") +
      scale_size(range = c(1.5, 6)) +
      scale_color_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
      geom_hline(yintercept = seq(4.5, 20.5, by = 4), linetype = "dashed", color = "black") +
      scale_y_discrete(labels = rep(c("Control", "B7-H4x4-1BB", "ICB", "Combi"), length(celltype_ann_level2))) +
      annotate("text",
               x     = n_genes + 0.8,
               y     = ct_pos,
               label = celltype_ann_level2,
               size  = 3, hjust = 0.5, fontface = "bold")

p <- VlnPlot(so.plot, group.by = 'ann_level2', split.by = 'hash.merged', features = genes, cols = hash.cols, pt.size = 0) & geom_boxplot(width = 0.1)
ggsave(p, filename = paste0(save_path, '/suppressive_vlnplot_LongAnnotation.png'), width = 10, height = 5)

Idents(so.cd4) <- so.cd4$ann_level2
markers <- FindMarkers(so.cd4, ident.1 = c('fTreg', 'fTreg.prolif.'), only.pos = T)

ggsave(paste0(save_path, '/suppressive_dotplot_LongAnnotation.png'), p2, width = 10, height = 5)
ggsave(paste0(save_path, '/suppressive_dotplot_LongAnnotation.pdf'), p2, width = 10, height = 5)
saveRDS(p2, paste0(save_path, '/suppressive_dotplot_LongAnnotation.rds'))
ggplot2pptx(p2, paste0(save_path, '/suppressive_dotplot_LongAnnotation.pptx'), width = 10, height = 5)

data <- p2$data
write_xlsx(data, paste0(save_path, "/excel_suppressive_dotplot_LongAnnotation_RawData.xlsx"))

############## 1.3.2 effector signature / ann_level_2
celltype_ann_level2 = c('CD4.Th1.like', 'fTreg', 'fTreg.prolif.', 'cTreg', 'cTreg.prolif.', 'Treg.LN.homing', 'Treg.IFN.response')
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
ggsave(p, file = paste0(save_path, '/supple7_LongAnnotation_signature_score_heatmap.png'), width = 10, height = 4)
ggsave(p, file = paste0(save_path, '/supple7_longAnnotation_signature_score_heatmap.pdf'), width = 10, height = 4)
ggplot2pptx(p, file = paste0(save_path, '/supple5_longAnnotation_signature_score_heatmap.pptx'), width = 10, height = 4)

raw_data_for_excel <- p$data
write_xlsx(raw_data_for_excel, paste0(save_path, "/excel_supple7_Heatmap_RawData_eff_exh.xlsx"))
