library(sceasy)
library(Seurat)
library(reticulate)
library(ggplot2)
library(magrittr)
library(data.table)
library(dplyr)
library(tidyr)
library(openxlsx)

source('/home/sjcho/yard/functions/R/save_ggplot2_to_ppt.R')
source('/home/sjcho/projects/utils.R')

save_path = '/home/sjcho/projects/4-1BB/publication_figures/new/publish_human_study/outs/main1_supple2'
set_name = 'CD8_ATLAS_healthy_solidTumor'

so <- readRDS('/home/sjcho/projects/4-1BB/public_data/CD8_atlas_analysis/outs/CD8_ATLAS_healthy_solidTumor_log_normalized.rds')
so.solid.tumor <- subset(so, cells = Cells(so)[so$meta_tissue_type == 'TIL' & so$treatment_status == 'Treatment-naive'])
so.solid.tumor$disease_control = ifelse(grepl('Control', so.solid.tumor$disease), 'Control', so.solid.tumor$disease %>% as.vector)

cold = c('Basal cell carinoma tumor', 'Breast cancer', 'Endometrial Carcinoma', 'Esophagus squamous cell carcinoma', 'Gallbladder carcinoma', 'low-grade gliomas', 'Ovarian Cancer', 'Pancreatic Cancer')
hot = c('Melanoma', 'HNSCC', 'Clear cell renal cell carcinoma', 'Gastric Cancer', 'Nasopharyngeal carcinoma', 'NSCLC', 'Squamous cell carcinoma tumor', 'Human thymic epithelial tumors (TET)')
not_classified = c('Metastatic colorectal cancer', 'Metastatic Cancer', 'Prostate Cancer', 'Cholangiocarcinoma')
# control = c('Control', 'Healthy', 'Healthy(Aging)', 'Healthy(CordBlood)', 'Young')
so.solid.tumor$disease_control = factor(so.solid.tumor$disease_control, 
                                        levels = c(cold, hot, not_classified)) # control
so.solid.tumor <- subset(so.solid.tumor, cells = Cells(so.solid.tumor)[!(so.solid.tumor$disease_control %in% c('Metastatic Cancer'))])

########### set order by TNFRSF9 mean in expression in disease_control
mean.exp.c <- c()
for (cancer_type in unique(so.solid.tumor$disease_control)) {
  mean.exp <- mean(expm1(so.solid.tumor@assays$RNA$data['TNFRSF9', Cells(so.solid.tumor)[so.solid.tumor$disease_control == cancer_type]]))
  mean.exp.c <- c(mean.exp.c, mean.exp)
}
names(mean.exp.c) <- unique(so.solid.tumor$disease_control)
mean.exp.c <- sort(mean.exp.c, decreasing = F)
so.solid.tumor$disease_control <- factor(as.vector(so.solid.tumor$disease_control), levels = names(mean.exp.c))

orange_to_red = c('#f9dacd', '#FCBBA1', '#FC9272', '#FB6A4A', '#EF3B2C', '#CB181D', '#A50F15', '#67000D')

p <- DotPlot(so.solid.tumor, features = c('TNFRSF9'), group.by = 'disease_control', scale = FALSE) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  scale_color_gradientn(colors = orange_to_red, limits = c(0, NA), breaks = c(0, 0.5, 1, 1.5))
ggsave(paste0(save_path, '/supple2_41BB_dotplot_by_disease_control_RdBu_ordered_noScale.png'), p, width = 6.5, height = 5)
ggsave(paste0(save_path, '/supple2_41BB_dotplot_by_disease_control_RdBu_ordered_noScale.pdf'), p, width = 6.5, height = 5)
ggplot2pptx(p, paste0(save_path, '/supple2_41BB_dotplot_by_disease_control_RdBu_ordered_noScale.pptx'), width = 6.5, height = 5)

p <- DotPlot(so.solid.tumor, features = c('TNFRSF9'), group.by = 'disease_control', scale = TRUE) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) + set_color_RdBu
ggsave(paste0(save_path, '/supple2_41BB_dotplot_by_disease_control_RdBu_ordered_Scale.png'), p, width = 6.5, height = 5)
ggsave(paste0(save_path, '/supple2_41BB_dotplot_by_disease_control_RdBu_ordered_Scale.pdf'), p, width = 6.5, height = 5)
ggplot2pptx(p, paste0(save_path, '/supple2_41BB_dotplot_by_disease_control_RdBu_ordered_Scale.pptx'), width = 6.5, height = 5)
saveRDS(p, file = paste0(save_path, '/supple2_41BB_dotplot_by_disease_control_RdBu_ordered_Scale.rds'))


p <- DotPlot(so.solid.tumor, features = c('TNFRSF9', 'ICOS', 'TNFRSF4', 'TNFRSF18'), group.by = 'disease_control', scale = FALSE) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        plot.title = element_text(size = 16, hjust = 0.5),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12)) +
  scale_color_gradientn(colors = orange_to_red)
ggsave(paste0(save_path, '/supple2_costimulatory_dotplot_by_disease_control_RdBu_ordered_noScale.png'), p, width = 8, height = 5)
ggsave(paste0(save_path, '/supple2_costimulatory_dotplot_by_disease_control_RdBu_ordered_noScale.pdf'), p, width = 8, height = 5)
ggplot2pptx(p, paste0(save_path, '/supple2_costimulatory_dotplot_by_disease_control_RdBu_ordered_noScale.pptx'), width = 8, height = 5)
saveRDS(p, file = paste0(save_path, '/supple2_costimulatory_dotplot_by_disease_control_RdBu_ordered_noScale.rds'))

wb <- createWorkbook()
addWorksheet(wb, "TNFRSF9")
writeData(wb, "TNFRSF9", p$data)
saveWorkbook(wb, file = paste0(save_path, '/supple2_whole_CD8.xlsx'), overwrite = TRUE)

####### supple2
cluster.groups <- list(
  CD8.naive = c(
    "Tn"
  ),
  CD8.Tcm = c(
    'Tcm', 'Early Tcm/Tem'
  ),
  CD8.Tem = c(
    'GZMK+ Tem', 'CMC1+ Temra', 'GNLY+ Temra'
  ),
  CD8.Trm = c(
    'CREM+ Trm', 'ITGB2+ Trm', 'ITGAE+ Trm'
  ),
  CD8.Teff = c(
    'ZNF683+ Teff'
  ),
  CD8.Tex = c(
    'XBP1+ Tex', 'ITGAE+ Tex', 'GZMK+ Tex'
  ),
  CD8.Tpex = c(
    'Tpex'
  ),
  CD8.ILTCK = c(
    'ILTCK'
  ),
  CD8.MAIT = c(
    'MAIT'
  ),
  CD8.MACF1_T = c(
    'MACF1+ T'
  ),
  CD8.Cycling_T = c(
    'Cycling T'
  )
)

so.solid.tumor$short_annotation <- 'not yet'
for (celltype in names(cluster.groups)) {
    so.solid.tumor@meta.data[Cells(so.solid.tumor)[so.solid.tumor$cell_subtype_3 %in% cluster.groups[[celltype]]], 'short_annotation'] = celltype
}
so.solid.tumor$short_annotation = factor(so.solid.tumor$short_annotation, levels = names(cluster.groups))

source('/home/sjcho/yard/functions/R/draw_single_gene_dotplot_by_two_metadata.R')
for (gene in c('TNFRSF9', 'ICOS', 'TNFRSF18', 'TNFRSF4')) { # 4-1BB, ICOS, GITR, OX40
  p <- draw_two_metadata_one_gene_expression_dotplot(so.solid.tumor, feature = gene, scale_expr = 'none',
                                                    row_metadata = 'disease_control', 
                                                    col_metadata = 'short_annotation')
  p <- p + theme(text = element_text(size = 20), 
                axis.text.x = element_text(angle = 90, hjust = 1))
  p <- p + scale_color_gradientn(colors = orange_to_red, limits = c(0, NA))
  ggsave(p, file = paste0(save_path, '/supple2_dotplot_', gene, '_disease_control_cell_subtype_3.png'), width = 7, height = 7)
  ggsave(p, file = paste0(save_path, '/supple2_dotplot_', gene, '_disease_control_cell_subtype_3.pdf'), width = 7, height = 7)
  ggplot2pptx(p, file = paste0(save_path, '/supple2_dotplot_', gene, '_disease_control_cell_subtype_3.pptx'), width = 7, height = 7)
}

wb <- createWorkbook()
for (gene in c('TNFRSF9', 'ICOS', 'TNFRSF18', 'TNFRSF4')) {
    p <- draw_two_metadata_one_gene_expression_dotplot_scaled(
        so.solid.tumor, feature = gene,
        row_metadata = 'disease_control',
        col_metadata = 'short_annotation'
    )
    addWorksheet(wb, gene)
    writeData(wb, gene, p$data)
}
saveWorkbook(wb, file = paste0(save_path, '/supple2_dotplot_BCDE_raw_data.xlsx'), overwrite = TRUE)