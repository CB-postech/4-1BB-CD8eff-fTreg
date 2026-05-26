### conda activate project_lung_exercise_R

library(Seurat)
library(magrittr)
library(data.table)
library(ggplot2)

save_path = '/home/sjcho/projects/4-1BB/publication_figures/new/publish_human_study/outs/main1_GSE156728'

source('/home/sjcho/yard/functions/R/draw_single_gene_dotplot_by_two_metadata.R')
source('/home/sjcho/projects/utils.R')
source('/home/sjcho/yard/functions/R/save_ggplot2_to_ppt.R')

meta.data <- fread('/home/sjcho/datas/public_data/GSE156728/metadata_and_rawfile/GSE156728_metadata.txt')
meta.data = data.frame(meta.data)
rownames(meta.data) = meta.data$cellID

base_dir <- "/home/sjcho/datas/public_data/GSE156728"
files <- list.files(base_dir, full.names = TRUE, recursive = FALSE)
files <- setdiff(files, c('/home/sjcho/datas/public_data/GSE156728/metadata_and_rawfile', '/home/sjcho/datas/public_data/GSE156728/data_from_zenodo'))

############### 1. load
all_counts <- list()
for (file in files) {   
    tissue <- sub("^GSE156728_([^_]+).*", "\\1", basename(file))
    CD4_or_CD8 <- sub(".*\\.(CD[48])\\.counts\\.txt$", "\\1", file)

    counts <- fread(file)
    counts <- as.data.frame(counts)
    rownames(counts) <- counts$V1
    counts$V1 <- NULL
    
    all_counts[[paste0(CD4_or_CD8, '_', tissue)]] <- counts
}

############## 2. Build seurat object
############## 2.1. find common genes
common_genes <- Reduce(intersect, lapply(all_counts, rownames))
print(common_genes)

############## 2.2. extract common genes
all_counts_common <- lapply(all_counts, function(x) x[common_genes, , drop = FALSE])
final_counts <- do.call(cbind, all_counts_common)
colnames(final_counts) <- unlist(lapply(all_counts_common, colnames))

dim(final_counts)
head(final_counts[, 1:5])

############## 2.3. extract common genes
so <- CreateSeuratObject(counts = final_counts, meta.data = meta.data[colnames(final_counts), ])

############## 3. Preprocessing
so <- NormalizeData(so)

############## 4. Visualization
############## 4.1. low resolution annotation
cluster.groups <- list(
  CD4.naive = c(
    "CD4.c01.Tn.TCF7",
    "CD4.c02.Tn.PASK",
    "CD4.c03.Tn.ADSL",
    "CD4.c04.Tn.il7r"
  ),
  CD4.memory = c(
    "CD4.c05.Tm.TNF",
    "CD4.c06.Tm.ANXA1",
    "CD4.c07.Tm.ANXA2",
    "CD4.c08.Tm.CREM",
    "CD4.c09.Tm.CCL5",
    "CD4.c10.Tm.CAPG",
    "CD4.c11.Tm.GZMA"
  ),
  CD4.Tem = c(
    "CD4.c12.Tem.GZMK",
    "CD4.c13.Temra.CX3CR1"
  ),
  CD4.Th17 = c(
    "CD4.c14.Th17.SLC4A10",
    "CD4.c15.Th17.IL23R"
  ),
  CD4.Tfh = c(
    "CD4.c16.Tfh.CXCR5"
  ),
  CD4.Th1 = c(
    "CD4.c17.TfhTh1.CXCL13"
  ),
  CD4.Treg = c(
    "CD4.c18.Treg.RTKN2",
    "CD4.c19.Treg.S1PR1",
    "CD4.c20.Treg.TNFRSF9",
    "CD4.c21.Treg.OAS1"
  ),
  CD4.ISG = c(
    "CD4.c22.ISG.IFIT1"
  ),
  CD4.MIX = c(
    "CD4.c23.Mix.NME1",
    "CD4.c24.Mix.NME2"
  ),
    CD8.naive = c(
    "CD8.c01.Tn.MAL"
  ),
  CD8.memory = c(
    "CD8.c02.Tm.IL7R",
    "CD8.c03.Tm.RPS12",
    "CD8.c04.Tm.CD52",
    "CD8.c17.Tm.NME1"
  ),
  CD8.Tem = c(
    "CD8.c05.Tem.CXCR5",
    "CD8.c06.Tem.GZMK",
    "CD8.c07.Temra.CX3CR1"
  ),
  CD8.NK.like = c(
    "CD8.c08.Tk.TYROBP",
    "CD8.c09.Tk.KIR2DL4"
  ),
  CD8.Trm = c(
    "CD8.c10.Trm.ZNF683"
  ),
  CD8.exhausted = c(
    "CD8.c11.Tex.PDCD1",
    "CD8.c12.Tex.CXCL13",
    "CD8.c13.Tex.myl12a",
    "CD8.c14.Tex.TCF7"
  ),
  CD8.ISG = c(
    "CD8.c15.ISG.IFIT1"
  ),
  MAIT = c(
    "CD8.c16.MAIT.SLC4A10"
  )
)

so$short_annotation <- 'not yet'
for (celltype in names(cluster.groups)) {
    so@meta.data[Cells(so)[so$meta.cluster %in% cluster.groups[[celltype]]], 'short_annotation'] = celltype
}
so$short_annotation = factor(so$short_annotation, levels = names(cluster.groups))

############## extract TIL
so.TIL <- subset(so, cells = Cells(so)[(so$loc == 'T') & !((so$cancerType %in% c('MM', 'BCL')))])
# without MM and BCL, because they are not solid tumors and may have different TME characteristics.

so.TIL$short_annotation <- 'not yet'
for (celltype in names(cluster.groups)) {
    so.TIL@meta.data[Cells(so.TIL)[so.TIL$meta.cluster %in% cluster.groups[[celltype]]], 'short_annotation'] = celltype
}
so.TIL$short_annotation = factor(so.TIL$short_annotation, levels = names(cluster.groups))
so.TIL$cancerType = factor(so.TIL$cancerType, levels = rev(c('ESCA', 'RC', 'UCEC', 'BC','PACA', 'THCA', 'OV', 'FTC')))

############## 1. 4-1BB across TIL and cancer types
orange_to_red = c('#f9dacd', '#FCBBA1', '#FC9272', '#FB6A4A', '#EF3B2C', '#CB181D', '#A50F15', '#67000D')

p1 <- DotPlot(so.TIL, features = 'TNFRSF9', group.by = 'short_annotation', scale = FALSE)
p1 <- p1 + scale_size(range = c(1, 7))
p1 <- p1 + scale_color_gradientn(
    colors = orange_to_red,
    limits = c(0, NA), 
    breaks = c(0, 0.5, 1, 1.5)
  )
ggsave(paste0(save_path, '/TNFRSF9_dotplot_by_short_annotation.png'), p1, width = 4.5, height = 5)
ggsave(paste0(save_path, '/TNFRSF9_dotplot_by_short_annotation.pdf'), p1, width = 4.5, height = 5)
ggplot2pptx(p1, paste0(save_path, '/TNFRSF9_dotplot_by_short_annotation.pptx'), width = 4.5, height = 5)

p2 <- DotPlot(so.TIL, features = 'TNFRSF9', group.by = 'cancerType', scale = TRUE)
p2 <- p2 + scale_size(range = c(1, 7))
p2 <- p2 + set_color_RdBu
ggsave(paste0(save_path, '/TNFRSF9_dotplot_by_cancerType.png'), p2, width = 4.5, height = 5)
ggsave(paste0(save_path, '/TNFRSF9_dotplot_by_cancerType.pdf'), p2, width = 4.5, height = 5)
ggplot2pptx(p2, paste0(save_path, '/TNFRSF9_dotplot_by_cancerType.pptx'), width = 4.5, height = 5)
saveRDS(p2, paste0(save_path, '/TNFRSF9_dotplot_by_cancerType.rds'))

cold_tumors = c('UCEC', 'PACA', 'OV', 'FTC', 'ESCA')
hot_tumors = c('RC')
Not_clear = c('THCA', 'BC')

so.TIL.exh = subset(so.TIL, subset = short_annotation %in% c('CD8.exhausted'))
mean.exp.c <- c()
for (cancer_type in as.vector(unique(so.TIL.exh$cancerType))) {
  mean.exp <- mean(expm1(so.TIL.exh@assays$RNA$data['TNFRSF9', Cells(so.TIL.exh)[so.TIL.exh$cancerType == cancer_type]]))
  mean.exp.c <- c(mean.exp.c, mean.exp)
}
names(mean.exp.c) <- as.vector(unique(so.TIL.exh$cancerType))
mean.exp.c <- sort(mean.exp.c, decreasing = F)

so.TIL.exh.Treg = subset(so.TIL, subset = short_annotation %in% c('CD4.Treg', 'CD8.exhausted'))
so.TIL.exh.Treg$short_annotation = factor(so.TIL.exh.Treg$short_annotation, levels = c('CD8.exhausted', 'CD4.Treg'))
so.TIL.exh.Treg$cancerType = factor(so.TIL$cancerType, levels = names(mean.exp.c))
p2 <- draw_two_metadata_one_gene_expression_dotplot(so.TIL.exh.Treg, 'TNFRSF9', 'cancerType', 'short_annotation', text_size = 5, scale_expr = 'column')
p2 <- p2 + set_color_RdBu
p2 <- p2 + coord_flip() + scale_size(range = c(1, 7))
p2 <- p2 + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(paste0(save_path, '/4-1BB_exh_Treg.png'), p2, width = 3.5, height = 5)
ggsave(paste0(save_path, '/4-1BB_exh_Treg.pdf'), p2, width = 3.5, height = 5)
ggplot2pptx(p2, paste0(save_path, '/4-1BB_exh_Treg.pptx'), width = 3.5, height = 5)

p <- combined_plot <- p1 + p2 + plot_layout(widths = c(1, 1.25))
ggsave(paste0(save_path, '/4-1BB_exh_Treg_combined.png'), p, width = 7, height = 7)
ggsave(paste0(save_path, '/4-1BB_exh_Treg_combined.pdf'), p, width = 7, height = 7)
ggplot2pptx(p, paste0(save_path, '/4-1BB_exh_Treg_combined.pptx'), width = 7, height = 7)
