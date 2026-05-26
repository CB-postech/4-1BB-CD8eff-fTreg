library(Seurat)
library(Matrix)
library(tidyr)
library(dplyr)
library(scRepertoire)
library(ggpubr)
library(patchwork)
library(broom)
library(tibble)
library(ggplotify)

source('/home/sjcho/projects/utils.R')
source('/home/sjcho/yard/functions/R/save_ggplot2_to_ppt.R')

save_path = '/home/sjcho/projects/4-1BB/publication_figures/new/publish_human_study/outs/supple8_fTreg_scVault'
data_dir <- "/home/sjcho/projects/4-1BB/public_data/scVault/Treg_relevance/outs/"

counts <- readMM("/home/sjcho/projects/4-1BB/public_data/scVault/outs/scVault_pre_post_datasets_tumor_R_NR_counts.mtx")

genes <- read.table(paste0(data_dir, "features.tsv"), header = FALSE, stringsAsFactors = FALSE)
barcodes <- read.table(paste0(data_dir, "barcodes.tsv"), header = FALSE, stringsAsFactors = FALSE)
metadata <- read.csv(paste0(data_dir, "metadata.csv"), row.names = 1) 
celltype <- read.csv('/home/sjcho/datas/public_data/scVault_curated_scRNA_scTCRseq/tumor_cluster_annotations.csv', row.names = 1)

colnames(counts) <- barcodes$V1
rownames(counts) <- genes$V1

so <- CreateSeuratObject(counts = counts)
so <- AddMetaData(object = so, metadata = metadata)
so$celltype <- celltype[Cells(so), 'cluster_annotations']

# CD8 T cells
CD8.T = c('C1_CD8_GZMK$^{+}$_1', 'C14_CD8_GZMB$^{+}$_2', 'C5_CD8_GZMK$^{+}$_2', 'C7_CD8_GZMB$^{+}$_1', 'C9_CD8_GZMB$^{+}$_Exhausted')

tcr_au <- read.csv("/home/sjcho/datas/public_data/scVault_curated_scRNA_scTCRseq/au_et_al/tcr_table_au_et_al.csv")
tcr_bassez_cohort_1 <- read.csv("/home/sjcho/datas/public_data/scVault_curated_scRNA_scTCRseq/bassez_et_al/cohort_1/tcr_table_bassez_et_al_cohort_1.csv")
tcr_bassez_cohort_2 <- read.csv("/home/sjcho/datas/public_data/scVault_curated_scRNA_scTCRseq/bassez_et_al/cohort_2/tcr_table_bassez_et_al_cohort_2.csv")
tcr_caushi <- read.csv("/home/sjcho/datas/public_data/scVault_curated_scRNA_scTCRseq/caushi_et_al/tcr_table_caushi_et_al.csv")
tcr_krishna <- read.csv("/home/sjcho/datas/public_data/scVault_curated_scRNA_scTCRseq/krishna_et_al/tcr_table_krishna_et_al.csv")
tcr_liu <- read.csv("/home/sjcho/datas/public_data/scVault_curated_scRNA_scTCRseq/liu_et_al/tcr_table_liu_et_al.csv")
tcr_luoma <- read.csv("/home/sjcho/datas/public_data/scVault_curated_scRNA_scTCRseq/luoma_et_al/tcr_table_luoma_et_al.csv")
tcr_pai <- read.csv("/home/sjcho/datas/public_data/scVault_curated_scRNA_scTCRseq/pai_et_al/tcr_table_pai_et_al.csv")
tcr_shiao <- read.csv("/home/sjcho/datas/public_data/scVault_curated_scRNA_scTCRseq/shiao_et_al/tcr_table_shiao_et_al.csv")
tcr_yost_BCC <- read.csv("/home/sjcho/datas/public_data/scVault_curated_scRNA_scTCRseq/yost_et_al/BCC/tcr_table_yost_et_al_BCC.csv")
tcr_yost_SCC <- read.csv("/home/sjcho/datas/public_data/scVault_curated_scRNA_scTCRseq/yost_et_al/SCC/tcr_table_yost_et_al_SCC.csv")
tcr_zhang <- read.csv("/home/sjcho/datas/public_data/scVault_curated_scRNA_scTCRseq/zhang_et_al/tcr_table_zhang_et_al.csv")

tcr_au$barcode = tcr_au$X; tcr_au$X <- NULL; tcr_au.filterd = tcr_au[tcr_au$barcode %in% Cells(so), ]
tcr_bassez_cohort_1$barcode = tcr_bassez_cohort_1$X; tcr_bassez_cohort_1$X <- NULL; tcr_bassez_cohort_1.filterd = tcr_bassez_cohort_1[tcr_bassez_cohort_1$barcode %in% Cells(so), ]
tcr_bassez_cohort_2$barcode = tcr_bassez_cohort_2$X; tcr_bassez_cohort_2$X <- NULL; tcr_bassez_cohort_2.filterd = tcr_bassez_cohort_2[tcr_bassez_cohort_2$barcode %in% Cells(so), ]
tcr_caushi$barcode = tcr_caushi$X; tcr_caushi$X <- NULL; tcr_caushi.filterd = tcr_caushi[tcr_caushi$barcode %in% Cells(so), ]
tcr_krishna$barcode = tcr_krishna$X; tcr_krishna$X <- NULL; tcr_krishna.filterd = tcr_krishna[tcr_krishna$barcode %in% Cells(so), ]
tcr_liu$barcode = tcr_liu$X; tcr_liu$X <- NULL; tcr_liu.filterd = tcr_liu[tcr_liu$barcode %in% Cells(so), ]
tcr_luoma$barcode = tcr_luoma$X; tcr_luoma$X <- NULL; tcr_luoma.filterd = tcr_luoma[tcr_luoma$barcode %in% Cells(so), ]
tcr_pai$barcode = tcr_pai$X; tcr_pai$X <- NULL; tcr_pai.filterd = tcr_pai[tcr_pai$barcode %in% Cells(so), ]
tcr_shiao$barcode = tcr_shiao$X; tcr_shiao$X <- NULL; tcr_shiao.filterd = tcr_shiao[tcr_shiao$barcode %in% Cells(so), ]
tcr_yost_BCC$barcode = tcr_yost_BCC$X; tcr_yost_BCC$X <- NULL; tcr_yost_BCC.filterd = tcr_yost_BCC[tcr_yost_BCC$barcode %in% Cells(so), ]
tcr_yost_SCC$barcode = tcr_yost_SCC$X; tcr_yost_SCC$X <- NULL; tcr_yost_SCC.filterd = tcr_yost_SCC[tcr_yost_SCC$barcode %in% Cells(so), ]
tcr_zhang$barcode = tcr_zhang$X; tcr_zhang$X <- NULL; tcr_zhang.filterd = tcr_zhang[tcr_zhang$barcode %in% Cells(so), ]

tcr_list = list('au' = tcr_au.filterd, 
     'bassez_1' = tcr_bassez_cohort_1.filterd,
     'bassez_2' = tcr_bassez_cohort_2.filterd,
     'caushi' = tcr_caushi.filterd,
     'krishna' = tcr_krishna.filterd,
     'liu' = tcr_liu.filterd,
     'luoma' = tcr_luoma.filterd,
     'pai' = tcr_pai.filterd,
     'shiao' = tcr_shiao.filterd,
     'yost_BCC' = tcr_yost_BCC.filterd,
     'yost_SCC' = tcr_yost_SCC.filterd,
     'zhang' = tcr_zhang.filterd)

tcr_merged <- rbind(tcr_au, tcr_bassez_cohort_1, tcr_bassez_cohort_2, tcr_caushi, tcr_krishna, tcr_liu, tcr_luoma, tcr_pai, tcr_shiao, tcr_yost_BCC, tcr_yost_SCC, tcr_zhang)

tcr_merged.filter <- tcr_merged[tcr_merged$barcode %in% Cells(so), ]
tcr_merged.filter.list <- list(tcr_merged.filter)
required_cols <- c("barcode", "chain", "contig_id", "v_gene", "d_gene", "j_gene", 
                   "c_gene", "cdr3", "cdr3_nt", "reads", "umis", "productive")

tcr_merged.filter.list <- lapply(tcr_merged.filter.list, function(df) {
  df %>%
    # 1. select required columns only
    select(any_of(required_cols)) %>%
    mutate(
      # 2. string to logical
      productive = case_when(
        tolower(as.character(productive)) == "true" ~ TRUE,
        tolower(as.character(productive)) == "false" ~ FALSE,
        TRUE ~ NA
      ),
      # 3. to numeric
      reads = as.numeric(reads),
      umis = as.numeric(umis),
      # 4. to character and trim whitespace
      chain = trimws(as.character(chain)),
      barcode = trimws(as.character(barcode))
    ) %>%
    as.data.frame()
})

combined.TCR <- combineTCR(list(tcr_merged.filter.list),
                           removeNA = FALSE, 
                           removeMulti = FALSE, 
                           filterMulti = TRUE)
so$TCR_clone_size <- 0
so$CTstrict <- 'Not Yet Assigned'
so$CTstrict[combined.TCR$S1$barcode] <- combined.TCR$S1$CTstrict

clone_size.strict <- table(combined.TCR$S1$CTstrict)
combined.TCR$S1$clone_size.strict <- as.integer(
  clone_size.strict[combined.TCR$S1$CTstrict]
)
so$TCR_clone_size[combined.TCR$S1$barcode] = combined.TCR$S1$clone_size.strict
so$is_expanded <- ifelse(so$TCR_clone_size > 1, 'Expanded', 'Not Expanded')
saveRDS(so, file = paste0(save_path, '/scVault_TIL_with_TCR.rds'))

########### Treg subset to find fTreg
so <- readRDS('/home/sjcho/projects/4-1BB/publication_figures/new/publish_human_study/outs/main6_CD8_41BB/scVault_TIL_with_TCR.rds')
so$celltype %>% table
so.Treg = subset(so, celltype %in% c('C2_Tregs'))
so.Treg <- log_normalize_harmony(so.Treg, save_path, 'Treg', pcs = 20, batch_name = 'study', nfeatures = 2000)

p <- DimPlot(so.Treg, group.by = 'seurat_clusters', label = TRUE)
ggsave(paste0(save_path, '/Treg_UMAP_by_cluster.png'), p, width = 6.5, height = 5)
p <- DimPlot(so.Treg, group.by = 'study', label = TRUE)
ggsave(paste0(save_path, '/Treg_UMAP_by_study.png'), p, width = 6.5, height = 5)

p <- DotPlot(so.Treg, features = c('FOXP3', 'IL2RA', 'IFNG', 'TNF', 'GZMB', 'GZMA', 'CCL4', 'CCL5', 'CD8A', 'CD4'), group.by = 'seurat_clusters', scale = TRUE) + set_color_RdBu
ggsave(paste0(save_path, '/Treg_dotplot_by_cluster_RdBu_ordered.png'), p, width = 6.5, height = 5)s
ggsave(paste0(save_path, '/Treg_dotplot_by_cluster_RdBu_ordered.pdf'), p, width = 6.5, height = 5)

### delete cluster 11, which is likely CD8 doublets
so.Treg.2 = subset(so.Treg, subset = seurat_clusters != 8)

so.Treg.2 <- log_normalize_harmony(so.Treg.2, save_path, 'Treg_2', pcs = 25, batch_name = 'study', nfeatures = 2000)
so.Treg.2 <- FindClusters(so.Treg.2, resolution = 1)

p <- DimPlot(so.Treg.2, group.by = 'seurat_clusters', label = FALSE, cols = c("#d58e2f", "#746fc9", "#54b84f", "#b95cc1", "#9cb93f", "darkred", 
"#63c499", "#cf444b", "#368868", "#ca5428", "#61a0d5", "#c4ad4e", "#bb6b8f", "#5c8537", "#cf8062", "#856f2e"))
p <- p + set_UMAP
ggsave(paste0(save_path, '/fig_fTreg2_UMAPcluster.png'), p, width = 5.5, height = 5)
saveRDS(as.ggplot(get_legend(p)), file = paste0(save_path, '/fig_fTreg2_UMAPcluster_legend.rds'))

p_legend <- get_legend(p)
ggplot2pptx(as_ggplot(p_legend), paste0(save_path, '/fig_Treg_UMAPcluster_legend.pptx'), width = 3, height = 4)

target_cells <- Cells(subset(so.Treg.2, subset = seurat_clusters == 5))
umap_df <- as.data.frame(Embeddings(so.Treg.2, "umap")[target_cells, ])
p <- DimPlot(so.Treg.2, cells.highlight = target_cells, cols.highlight = 'darkred', sizes.highlight = 0.05, cols = 'lightgray') +
  geom_density_2d(data = umap_df, aes(x = umap_1, y = umap_2), color = "black", linewidth = 0.5, bins = 7)
p <- p + NoLegend() + set_UMAP
ggsave(paste0(save_path, '/fig_Treg_UMAPcluster_highlight_cluster5.png'), p, width = 5, height = 5)

p <- FeaturePlot(so.Treg.2, features = c('CD4', 'FOXP3', 'CD8A', 'GNLY', 'GZMB', 'NKG7'), ncol = 3, order = T) & set_UMAP_featureplot
p <- p & scale_color_gradientn(colors = brewer.pal(9, "Reds"))
# decraese the size of legends
p <- p & theme(legend.key.size = unit(0.5, "cm"), legend.text = element_text(size = 6), legend.title = element_text(size = 8))
saveRDS(p, file = paste0(save_path, '/fig_Treg_2_featureplot.rds'))

for (genes in c('CD4', 'FOXP3', 'CD8A', 'GNLY', 'GZMB', 'NKG7')) {
  p <- FeaturePlot(so.Treg.2, features = genes, ncol = 1, order = T) + set_UMAP_featureplot +
    scale_color_gradientn(colors = brewer.pal(9, "Reds")) +
    theme(legend.key.size = unit(0.5, "cm"), legend.text = element_text(size = 6), legend.title = element_text(size = 8))
  p_leg <- p + theme_void() + theme(legend.position = "right")
  saveRDS(p_leg, file = paste0(save_path, '/fig_Treg_2_featureplot_', genes, '_legend.rds'))
}

ggsave(paste0(save_path, '/fig_Treg_2_featureplot.png'), p, width = 10, height = 6)

so.fTreg <- subset(so.Treg.2, subset = seurat_clusters == 5)
saveRDS(so.fTreg, file = paste0(save_path, '/so.fTreg.rds'))

so.fTreg <- readRDS('/home/sjcho/projects/4-1BB/publication_figures/new/publish_human_study/outs/supple8_fTreg_scVault/so.fTreg.rds')
so.wo.doublets <- subset(so, cells = c(setdiff(Cells(so), Cells(so.Treg)), Cells(so.Treg.2)))
so.wo.doublets$celltype_fTreg <- ifelse(Cells(so.wo.doublets) %in% Cells(so.fTreg), 'fTreg', so.wo.doublets$celltype)

################### cytotoxicity & 4-1BB response level / ssGSEA
library(msigdbr)
# BioCarta
biocarta_41bb <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:BIOCARTA") |>
  dplyr::filter(gs_name == "BIOCARTA_41BB_PATHWAY")
# GO:BP
gobp_cytotox <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") |>
  dplyr::filter(gs_name == "GOBP_LEUKOCYTE_MEDIATED_CYTOTOXICITY")

gobp_Tactivation <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") |>
  dplyr::filter(gs_name == "GOBP_T_CELL_ACTIVATION")

GOCC_CYTOLYTIC_GRANULE <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:CC") |>
  dplyr::filter(gs_name == "GOCC_CYTOLYTIC_GRANULE")

ctl_effector <- c(
  "PRF1", "GZMA", "GZMB", "GZMH", "GZMK", "GZMM", "GNLY",
  "NKG7", "FASLG", "TNF"
)

# degranulation & trafficking
# RAB27A, STX11, STX7, LAMP1

# GOCC_CYTOLYTIC_GRANULE

intersect(biocarta_41bb$gene_symbol, gobp_cytotox$gene_symbol)
intersect(biocarta_41bb$gene_symbol, ctl_effector)
intersect(biocarta_41bb$gene_symbol, GOCC_CYTOLYTIC_GRANULE$gene_symbol)

######## sample level GSVA/ssGSEA scores by cell type
save_path = '/home/sjcho/projects/4-1BB/publication_figures/new/publish_human_study/outs/supple8_fTreg_scVault/CPM_normalized_ssGSEA_41BB_custom_cytotoxic_test_original_anno'

so.fTreg <- readRDS('/home/sjcho/projects/4-1BB/publication_figures/new/publish_human_study/outs/supple8_fTreg_scVault/so.fTreg.rds')
so.wo.doublets <- subset(so, cells = c(setdiff(Cells(so), Cells(so.Treg)), Cells(so.Treg.2)))
so.wo.doublets$celltype_fTreg <- ifelse(Cells(so.wo.doublets) %in% Cells(so.fTreg), 'fTreg', so.wo.doublets$celltype)
so.wo.doublets$hot_cold = ifelse(so.wo.doublets$cancer_type %in% c('HNSCC', 'RCC', 'NSCLC', 'BCC'), 'hot', 'cold')

library(GSVA)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# gene sets (from your code)
gene_sets <- list(
  `41BB_pathway` = biocarta_41bb$gene_symbol,
  cytotoxic = ctl_effector
)

# 1. Pseudobulk 생성
make_pseudobulk <- function(so, n_filter) {
  cell_counts <- so@meta.data |>
    dplyr::count(sample, celltype_fTreg) |>
    dplyr::filter(n >= n_filter)
  
  cells_keep <- so@meta.data |>
    tibble::rownames_to_column("barcode") |>
    dplyr::inner_join(cell_counts, by = c("sample", "celltype_fTreg")) |>
    dplyr::pull(barcode)
  
  if (length(cells_keep) == 0) return(NULL)
  
  bulk <- subset(so, cells = cells_keep) |>
    AggregateExpression(group.by = c("sample", "celltype_fTreg"), return.seurat = TRUE)
  NormalizeData(bulk, normalization.method = "LogNormalize", scale.factor = 1e6)
}

# 2. ssGSEA scoring
run_ssgsea <- function(bulk, gene_sets) {
  mat <- GetAssayData(bulk, slot = "data") |> as.matrix()
  scores <- gsva(mat, gene_sets, method = "ssgsea") |> t() |> as.data.frame()
  scores$id <- rownames(scores)
  scores$celltype_fTreg <- sub(".*_", "", scores$id)
  scores$sample <- sub("_[^_]+$", "", scores$id)
  scores$method <- "ssgsea"
  scores
}

# 3. Boxplot
plot_boxplot <- function(df_scores, save_path, tag, n_filter) {
  df_long <- df_scores |>
    pivot_longer(cols = c(`41BB_pathway`, cytotoxic), names_to = "gene_set", values_to = "score")
  
  p <- ggplot(df_long, aes(x = celltype_fTreg, y = score, fill = celltype_fTreg)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.2, size = 0.1, alpha = 1) +
    facet_wrap(~ gene_set, scales = "free_y") +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          legend.position = "none") +
    labs(x = NULL, y = "SSGSEA Score", title = tag)
  
  ggsave(p, filename = paste0(save_path, "/boxplot_celltype_41BB_cytotoxic_ssgsea_", tag, "_min", n_filter, ".png"),
         width = 12, height = 6)
}

# 4. Metadata 붙이기
attach_meta <- function(df_scores, so) {
  meta_map <- so@meta.data %>% distinct(sample, treatment, response)
  meta_map$sample_agg <- gsub("_", "-", meta_map$sample)
  
  df <- df_scores %>% left_join(meta_map, by = c("sample" = "sample_agg"))
  df$response <- factor(df$response, levels = c("NR", "R"))
  df$treatment <- factor(df$treatment, levels = c("Pre", "Post"))
  df
}

# 5. Scatter plot (per celltype)
plot_scatter <- function(df_scores, save_path, tag, n_filter) {
  for (ct in unique(df_scores$celltype_fTreg)) {
    print(paste("Processing:", tag, ct))
    df_ct <- df_scores %>% filter(celltype_fTreg == ct)
    
    cor_labels <- df_ct %>%
      group_by(treatment, response) %>%
      summarise(
        n = sum(is.finite(`41BB_pathway`) & is.finite(cytotoxic)),
        ct = ifelse(n >= 3, list(cor.test(`41BB_pathway`, cytotoxic)), list(NULL)),
        .groups = "drop"
      ) %>%
      mutate(label = sapply(ct, \(x) if (is.null(x)) NA_character_ else sprintf("r=%.3f\np=%.2e", x$estimate, x$p.value)))
    
    df_smooth <- df_ct %>%
      semi_join(cor_labels %>% filter(n >= 3), by = c("treatment", "response"))
    
    xl <- range(df_ct$`41BB_pathway`, na.rm = TRUE)
    yl <- range(df_ct$cytotoxic, na.rm = TRUE)
    
    p <- ggplot(df_ct, aes(x = `41BB_pathway`, y = cytotoxic)) +
      geom_point(size = 2, alpha = 0.7) +
      geom_smooth(data = df_smooth, method = "lm", se = TRUE, color = "black") +
      geom_text(data = cor_labels %>% filter(!is.na(label)), aes(label = label),
                x = -Inf, y = Inf, hjust = -0.1, vjust = 1.3, size = 3.5) +
      facet_grid(treatment ~ response) +
      xlim(xl) + ylim(yl) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
      labs(title = paste0(ct, " (", tag, ")"),
           x = "SSGSEA 41BB pathway",
           y = "SSGSEA Cytotoxic")
    
    ggsave(p, filename = paste0(save_path, "/scatter_41BB_vs_cytotoxic_trmt_resp_ssgsea_", ct, "_", tag, "_min", n_filter, ".png"), width = 10, height = 8)
    saveRDS(p, file = paste0(save_path, "/scatter_41BB_vs_cytotoxic_trmt_resp_ssgsea_", ct, "_", tag, "_min", n_filter, ".rds"))
  }
}

save_path = '/home/sjcho/projects/4-1BB/publication_figures/new/publish_human_study/outs/supple8_fTreg_scVault/CPM_normalized_ssGSEA_41BB_custom_cytotoxic_test'
so.wo.doublets$hot_cold <- ifelse(so.wo.doublets$cancer_type %in% c('HNSCC', 'RCC', 'NSCLC', 'BCC'), 'hot', 'cold')

for (n_filter in c(25, 27, 29, 31, 33, 35)) {
  for (hc in c("hot", "cold")) {
    so_sub <- subset(so.wo.doublets, hot_cold == hc)
    
    bulk <- make_pseudobulk(so_sub, n_filter)
    if (is.null(bulk)) next
    
    df_scores <- run_ssgsea(bulk, gene_sets)
    
    plot_boxplot(df_scores, save_path, hc, n_filter)
    
    df_scores <- attach_meta(df_scores, so_sub)
    
    plot_scatter(df_scores, save_path, hc, n_filter)
  }
}

########################## No R, NR, Pre, Post
save_path = '/home/sjcho/projects/4-1BB/publication_figures/new/publish_human_study/outs/supple8_fTreg_scVault/CPM_normalized_ssGSEA_41BB_custom_cytotoxic_test_No_R_NR_Pre_Post'
plot_scatter_simple <- function(df_scores, save_path, tag, n_filter) {
  for (ct in unique(df_scores$celltype_fTreg)) {
    print(paste("Processing:", tag, ct))
    df_ct <- df_scores %>% filter(celltype_fTreg == ct)
    
    n_finite <- sum(is.finite(df_ct$`41BB_pathway`) & is.finite(df_ct$cytotoxic))
    
    if (n_finite >= 3) {
      ct_res <- cor.test(df_ct$`41BB_pathway`, df_ct$cytotoxic)
      label <- sprintf("r=%.3f\np=%.2e", ct_res$estimate, ct_res$p.value)
    } else {
      label <- NULL
    }
    
    xl <- range(df_ct$`41BB_pathway`, na.rm = TRUE)
    yl <- range(df_ct$cytotoxic, na.rm = TRUE)
    
    p <- ggplot(df_ct, aes(x = `41BB_pathway`, y = cytotoxic)) +
      geom_point(size = 2, alpha = 0.7) +
      xlim(xl) + ylim(yl) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
      labs(title = paste0(ct, " (", tag, ")"),
           x = "SSGSEA 41BB pathway",
           y = "SSGSEA Cytotoxic")
    
    if (n_finite >= 3) {
      p <- p +
        geom_smooth(method = "lm", se = TRUE, color = "black") +
        annotate("text", x = -Inf, y = Inf, label = label,
                 hjust = -0.1, vjust = 1.3, size = 3.5)
    }
    
    ggsave(p, filename = paste0(save_path, "/scatter_41BB_vs_cytotoxic_ssgsea_", ct, "_", tag, "_min", n_filter, ".png"), width = 6, height = 5)
    saveRDS(p, file = paste0(save_path, "/scatter_41BB_vs_cytotoxic_ssgsea_", ct, "_", tag, "_min", n_filter, ".rds"))
  }
}

n_filter = 30
for (hc in c("hot", "cold")) {
  so_sub <- subset(so.wo.doublets, hot_cold == hc)
  bulk <- make_pseudobulk(so_sub, n_filter)
  if (is.null(bulk)) next
  df_scores <- run_ssgsea(bulk, gene_sets)
  plot_boxplot(df_scores, save_path, hc, n_filter)
  plot_scatter_simple(df_scores, save_path, hc, n_filter)
}

collect_cor <- function(df_scores, tag) {
  results <- list()
  for (ct in unique(df_scores$celltype_fTreg)) {
    df_ct <- df_scores %>% filter(celltype_fTreg == ct)
    n_finite <- sum(is.finite(df_ct$`41BB_pathway`) & is.finite(df_ct$cytotoxic))
    if (n_finite >= 3) {
      ct_res <- cor.test(df_ct$`41BB_pathway`, df_ct$cytotoxic)
      results[[ct]] <- data.frame(celltype = ct, group = tag,
                                  r = ct_res$estimate, p = ct_res$p.value)
    }
  }
  bind_rows(results)
}

# 2. 수집
cor_all <- list()
for (hc in c("hot", "cold")) {
  so_sub <- subset(so.wo.doublets, hot_cold == hc)
  bulk <- make_pseudobulk(so_sub, n_filter)
  if (is.null(bulk)) next
  df_scores <- run_ssgsea(bulk, gene_sets)
  cor_all[[hc]] <- collect_cor(df_scores, hc)
}
cor_df <- bind_rows(cor_all)
cor_df$neglog10p <- -log10(cor_df$p)
cor_df$group <- factor(cor_df$group, levels = c("hot", "cold"))

cor_df$sig <- ifelse(cor_df$p < 0.05, "*", "")

p <- ggplot(cor_df, aes(x = celltype, y = group, size = neglog10p, color = r)) +
  geom_point() +
  geom_text(aes(label = sig), color = "black", size = 5, vjust = 0.35, show.legend = FALSE) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                        name = "Pearson r") +
  scale_size_continuous(name = expression(-log[10](p))) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = NULL, y = NULL)

ggsave(p, filename = paste0(save_path, "/dotplot_cor_41BB_cytotoxic_hot_cold_min", n_filter, ".png"),
      width = 10, height = 7)
