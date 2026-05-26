### ab normalization based on https://satijalab.org/seurat/articles/multimodal_vignette.html
### HTOdemux based on https://satijalab.org/seurat/articles/hashing_vignette.html
### conda activate project_lung_exercise_R

library(Seurat)
library(SingleCellExperiment)
library(magrittr)
library(ggplot2)
library(reshape2)

save_path = '/home/sjcho/projects/4-1BB/DA240060/outs/QC/3.demultiplexing'

# sce.p1 <- readRDS('/home/sjcho/projects/4-1BB/DA240060/outs/QC/2.lowQuliatyQC/filtered_pool1_mt15.rds')
# sce.p2 <- readRDS('/home/sjcho/projects/4-1BB/DA240060/outs/QC/2.lowQuliatyQC/filtered_pool2_lower600_umi1000_mt15.rds')
sce.p1 <- readRDS('/home/sjcho/projects/4-1BB/DA240060/outs/QC/2.lowQuliatyQC/filtered_pool1.rds')
sce.p2 <- readRDS('/home/sjcho/projects/4-1BB/DA240060/outs/QC/2.lowQuliatyQC/filtered_pool2_lower600_umi2000_mt5.rds')

########## pool 1
count.p1 <- Read10X('/home/sjcho/projects/4-1BB/DA240060/DA240060_41BB_20250204_pool1/outs/multi/count/raw_feature_bc_matrix')
count.RNA.p1 <- count.p1$'Gene Expression'[, colnames(sce.p1)]; colnames(count.RNA.p1) <- paste0('pool1_', colnames(count.RNA.p1))
count.Ab.p1 <- count.p1$'Antibody Capture'[, colnames(sce.p1)]; colnames(count.Ab.p1) <- paste0('pool1_', colnames(count.Ab.p1))

so.p1 <- CreateSeuratObject(counts = count.RNA.p1, project = "pool1", min.cells = 0, min.features = 0)
hashtag.p1 <- CreateAssay5Object(counts = count.Ab.p1[c('control_T_totalseq_C301', 'bsAb_T_totalseq_C302', 'ICB_T_totalseq_C303'), ] )
cite.p1 <- CreateAssay5Object(counts = count.Ab.p1[c('PD1_totalseq_C_0004', 'Tim3_totalseq_C_0003', '41BB_totalseq_C_0194'), ])

# abc stand for AntiBody Capture
so.p1[['hashtag']] <- hashtag.p1
so.p1[['cite']] <- cite.p1

DefaultAssay(so.p1) <- 'RNA'
so.p1 <- NormalizeData(so.p1)
so.p1 <- FindVariableFeatures(so.p1)
so.p1 <- ScaleData(so.p1, features = VariableFeatures(so.p1))

DefaultAssay(so.p1) <- 'hashtag'
so.p1 <- NormalizeData(so.p1, normalization.method = "CLR", margin = 1)

### Hashtag histogram
data <- so.p1@assays$hashtag$data
fill.list = list(
    'control-T-totalseq-C301' = 'darkblue',
    'bsAb-T-totalseq-C302' = 'darkred',
    'ICB-T-totalseq-C303' = 'forestgreen'
)
cutoff.list = list(
    'control-T-totalseq-C301' = 0.8,
    'bsAb-T-totalseq-C302' = 1,
    'ICB-T-totalseq-C303' = 0.8
)

for (ab in c('control-T-totalseq-C301', 'bsAb-T-totalseq-C302', 'ICB-T-totalseq-C303')) {
  df <- data.frame(Value = as.numeric(data[ab, ]))
  sample_name <- rownames(data)[ab]
  
  p <- ggplot(df, aes(x = Value)) +
    geom_histogram(binwidth = diff(range(df$Value)) / 100,
                   fill    = fill.list[[ab]],
                   color   = "white") +
    ggtitle(ab) +
    geom_vline(xintercept = cutoff.list[[ab]], color = "black", linetype = 'dotted', size = 1) +
    theme_classic()
  ggsave(p, file = paste0(save_path, '/pool1_histogram_', ab, '.png'), width = 8, height = 6)
}

### HTO demultiplexing
so.p1 <- HTODemux(so.p1, assay = "hashtag", positive.quantile = 0.99)
table(so.p1$hashtag_classification.global)

Idents(so.p1) <- "hashtag_classification.global"
p <- VlnPlot(so.p1, features = "nCount_RNA", pt.size = 0, log = TRUE, cols = c('gray30', 'gray30', 'gray30')) + NoLegend() + geom_boxplot(fill = 'white', width = 0.2)
ggsave(p, file = paste0(save_path, '/HTODemux_nCountRNA_poo1.png'), width = 9, height = 6)

saveRDS(so.p1, file = paste0(save_path, '/pool1.rds'))
# then use hasd.ID

## raw count, whole droplets
data.full <- count.p1$'Antibody Capture'[c('control_T_totalseq_C301', 'bsAb_T_totalseq_C302', 'ICB_T_totalseq_C303'), ]
fill.list = list(
    'control_T_totalseq_C301' = 'darkblue',
    'bsAb_T_totalseq_C302' = 'darkred',
    'ICB_T_totalseq_C303' = 'forestgreen'
)

for (ab in c('control_T_totalseq_C301', 'bsAb_T_totalseq_C302', 'ICB_T_totalseq_C303')) {
  df <- data.frame(Value = log10(as.numeric(data.full[ab, ]) + 1))
  
  p <- ggplot(df, aes(x = Value)) +
    geom_histogram(binwidth = diff(range(df$Value)) / 250,
                   fill    = fill.list[[ab]],
                   color   = "white") +
    ylim(0, 10000) +
    ggtitle(ab) +
#     geom_vline(xintercept = cutoff.list[[ab]], color = "black", linetype = 'dotted', size = 1) +
    theme_classic()
  ggsave(p, file = paste0(save_path, '/pool1_histogram_', ab, '_whole_cell_count.png'), width = 8, height = 6)
}

# so.p1$demultipexing <- 'Not Assigned'
# # triplet
# so.p1$demultipexing[so.p1$abc$data['control-T-totalseq-C301', ] > cutoff.list[['control-T-totalseq-C301']] & 
#                 so.p1$abc$data['bsAb-T-totalseq-C302', ] > cutoff.list[['bsAb-T-totalseq-C302']] &
#                 so.p1$abc$data['ICB-T-totalseq-C303', ] > cutoff.list[['ICB-T-totalseq-C303']]] <- 'Triplet_control_bsAb_ICB'

# # doublet
# so.p1$demultipexing[so.p1$abc$data['control-T-totalseq-C301', ] > cutoff.list[['control-T-totalseq-C301']] & 
#                 so.p1$abc$data['bsAb-T-totalseq-C302', ] > cutoff.list[['bsAb-T-totalseq-C302']] &
#                 so.p1$abc$data['ICB-T-totalseq-C303', ] <= cutoff.list[['ICB-T-totalseq-C303']]] <- 'Doublet_control_bsAb'
# so.p1$demultipexing[so.p1$abc$data['control-T-totalseq-C301', ] > cutoff.list[['control-T-totalseq-C301']] &
#                 so.p1$abc$data['bsAb-T-totalseq-C302', ] <= cutoff.list[['bsAb-T-totalseq-C302']] &
#                 so.p1$abc$data['ICB-T-totalseq-C303', ] > cutoff.list[['ICB-T-totalseq-C303']]] <- 'Doublet_control_ICB'
# so.p1$demultipexing[so.p1$abc$data['control-T-totalseq-C301', ] <= cutoff.list[['control-T-totalseq-C301']] &
#                 so.p1$abc$data['bsAb-T-totalseq-C302', ] > cutoff.list[['bsAb-T-totalseq-C302']] &
#                 so.p1$abc$data['ICB-T-totalseq-C303', ] > cutoff.list[['ICB-T-totalseq-C303']]] <- 'Doublet_bsAb_ICB'

# # singlet
# so.p1$demultipexing[so.p1$abc$data['control-T-totalseq-C301', ] > cutoff.list[['control-T-totalseq-C301']] & 
#                 so.p1$abc$data['bsAb-T-totalseq-C302', ] <= cutoff.list[['bsAb-T-totalseq-C302']] &
#                 so.p1$abc$data['ICB-T-totalseq-C303', ] <= cutoff.list[['ICB-T-totalseq-C303']]] <- 'Singlet_control'
# so.p1$demultipexing[so.p1$abc$data['control-T-totalseq-C301', ] <= cutoff.list[['control-T-totalseq-C301']] &
#                 so.p1$abc$data['bsAb-T-totalseq-C302', ] > cutoff.list[['bsAb-T-totalseq-C302']] &
#                 so.p1$abc$data['ICB-T-totalseq-C303', ] <= cutoff.list[['ICB-T-totalseq-C303']]] <- 'Singlet_bsAb'
# so.p1$demultipexing[so.p1$abc$data['control-T-totalseq-C301', ] <= cutoff.list[['control-T-totalseq-C301']] &
#                 so.p1$abc$data['bsAb-T-totalseq-C302', ] <= cutoff.list[['bsAb-T-totalseq-C302']] &
#                 so.p1$abc$data['ICB-T-totalseq-C303', ] > cutoff.list[['ICB-T-totalseq-C303']]] <- 'Singlet_ICB'

########## pool 2
count.p2 <- Read10X('/home/sjcho/projects/4-1BB/DA240060/DA240060_41BB_20250204_pool2/outs/multi/count/raw_feature_bc_matrix')
count.RNA.p2 <- count.p2$'Gene Expression'[, colnames(sce.p2)]; colnames(count.RNA.p2) <- paste0('pool2_', colnames(count.RNA.p2))
count.Ab.p2 <- count.p2$'Antibody Capture'[, colnames(sce.p2)]; colnames(count.Ab.p2) <- paste0('pool2_', colnames(count.Ab.p2))

so.p2 <- CreateSeuratObject(counts = count.RNA.p2, project = "pool2", min.cells = 0, min.features = 0)
hashtag.p2 <- CreateAssay5Object(counts = count.Ab.p2[c('Combi_T_totalseq_C304', 'MC38_T_totalseq_C305'), ] )
cite.p2 <- CreateAssay5Object(counts = count.Ab.p2[c('PD1_totalseq_C_0004', 'Tim3_totalseq_C_0003', '41BB_totalseq_C_0194'), ])

so.p2[['hashtag']] <- hashtag.p2
so.p2[['cite']] <- cite.p2

DefaultAssay(so.p2) <- 'RNA'
so.p2 <- NormalizeData(so.p2)
so.p2 <- FindVariableFeatures(so.p2)
so.p2 <- ScaleData(so.p2, features = VariableFeatures(so.p2))

DefaultAssay(so.p2) <- 'hashtag'
so.p2 <- NormalizeData(so.p2, normalization.method = "CLR", margin = 1)

### Hashtag histogram
### raw count
data <- so.p2@assays$hashtag$count
data <- log10(data + 1)
fill.list = list(
    'Combi-T-totalseq-C304' = 'darkblue',
    'MC38-T-totalseq-C305' = 'darkred'
)
cutoff.list = list(
    'Combi-T-totalseq-C304' = 2.4,
    'MC38-T-totalseq-C305' = 2.75
)
for (ab in c('Combi-T-totalseq-C304', 'MC38-T-totalseq-C305')) {
  df <- data.frame(Value = as.numeric(data[ab, ]))
  sample_name <- rownames(data)[ab]
  
  p <- ggplot(df, aes(x = Value)) +
    geom_histogram(binwidth = diff(range(df$Value)) / 250,
                   fill    = fill.list[[ab]],
                   color   = "white") +
    xlim(1.5, 3.5) +
    ggtitle(ab) +
    geom_vline(xintercept = cutoff.list[[ab]], color = "black", linetype = 'dotted', size = 1) +
    theme_classic()
  ggsave(p, file = paste0(save_path, '/pool2_histogram_lower600_umi2000_', ab, '_count.png'), width = 8, height = 6)
}

df <- data.frame(
  x = as.numeric(data["Combi-T-totalseq-C304", ]),
  y = as.numeric(data["MC38-T-totalseq-C305", ]),
  sample = colnames(data)
)
p <- ggplot(df, aes(x = x, y = y)) +
  geom_point(alpha = 0.7, size = 3) +
  theme_classic() +
  labs(x = "Combi-T-totalseq-C304", y = "MC38-T-totalseq-C305")
ggsave(p, file = paste0(save_path, '/pool2_logCount_2D_lower600_umi2000_.png'), width = 7, height = 6)

# colnames(data)[data['Combi-T-totalseq-C304', ] < 2.4 & data['MC38-T-totalseq-C305', ] < 2.8] -> not.assigned
colnames(data)[data['Combi-T-totalseq-C304', ] > 2.4] -> combi.cells
colnames(data)[data['Combi-T-totalseq-C304', ] < 2.4] -> mc38.cells
# colnames(data)[data['Combi-T-totalseq-C304', ] > 2.4 & data['MC38-T-totalseq-C305', ] > 2.8] -> doublet

# colnames(data)[data['Combi-T-totalseq-C304', ] < 2.4 & data['MC38-T-totalseq-C305', ] < 2.8] -> not.assigned
# colnames(data)[data['Combi-T-totalseq-C304', ] > 2.4 & data['MC38-T-totalseq-C305', ] < 2.8] -> combi.cells
# colnames(data)[data['Combi-T-totalseq-C304', ] < 2.4 & data['MC38-T-totalseq-C305', ] > 2.8] -> mc38.cells
# colnames(data)[data['Combi-T-totalseq-C304', ] > 2.4 & data['MC38-T-totalseq-C305', ] > 2.8] -> doublet

so.p2$hash.manual <- 'negative'
so.p2$hash.manual[combi.cells] <- 'Combi-T-totalseq-C304'
so.p2$hash.manual[mc38.cells] <- 'MC38-T-totalseq-C305'
# so.p2$hash.manual[doublet] <- 'doublet'

# ### normalized count
# data <- so.p2@assays$hashtag$data
# fill.list = list(
#     'Combi-T-totalseq-C304' = 'darkblue',
#     'MC38-T-totalseq-C305' = 'darkred'
# )
# cutoff.list = list(
#     'Combi-T-totalseq-C304' = 0.4,
#     'MC38-T-totalseq-C305' = 0.8
# )
# for (ab in c('Combi-T-totalseq-C304', 'MC38-T-totalseq-C305')) {
#   df <- data.frame(Value = as.numeric(data[ab, ]))
#   sample_name <- rownames(data)[ab]
  
#   p <- ggplot(df, aes(x = Value)) +
#     geom_histogram(binwidth = diff(range(df$Value)) / 250,
#                    fill    = fill.list[[ab]],
#                    color   = "white") +
#     ggtitle(ab) +
#     xlim(0, 2) + 
#     geom_vline(xintercept = cutoff.list[[ab]], color = "black", linetype = 'dotted', size = 1) +
#     theme_classic()
#   ggsave(p, file = paste0(save_path, '/pool2_histogram_lower600_', ab, '_normalized.png'), width = 8, height = 6)
# }

# colnames(data)[data['Combi-T-totalseq-C304', ] < 0.4 & data['MC38-T-totalseq-C305', ] < 0.8] -> not.assigned
# colnames(data)[data['Combi-T-totalseq-C304', ] > 0.4 & data['MC38-T-totalseq-C305', ] < 0.8] -> combi.cells
# colnames(data)[data['Combi-T-totalseq-C304', ] < 0.4 & data['MC38-T-totalseq-C305', ] > 0.8] -> mc38.cells
# colnames(data)[data['Combi-T-totalseq-C304', ] > 0.4 & data['MC38-T-totalseq-C305', ] > 0.8] -> doublet

# so.p2$hash.manual.norm <- 'negative'
# so.p2$hash.manual.norm[combi.cells] <- 'Combi-T-totalseq-C304'
# so.p2$hash.manual.norm[mc38.cells] <- 'MC38-T-totalseq-C305'
# so.p2$hash.manual.norm[doublet] <- 'doublet'

# df <- data.frame(
#   x = as.numeric(data["Combi-T-totalseq-C304", ]),
#   y = as.numeric(data["MC38-T-totalseq-C305", ]),
#   sample = colnames(data)
# )
# p <- ggplot(df, aes(x = x, y = y)) +
#   geom_point(alpha = 0.7, size = 3) +
#   theme_classic() +
#   labs(x = "Combi-T-totalseq-C304", y = "MC38-T-totalseq-C305")
# ggsave(p, file = paste0(save_path, '/pool2_2D_normalized_lower500.png'), width = 7, height = 6)

### for yjjeong demultiplexing
write.csv(t(so.p2@assays$hashtag$data), file = paste0(save_path, '/pool2_lower600_UMI2000_normalized_Hashtag.csv'))
write.csv(t(log10(so.p2@assays$hashtag$count + 1)), file = paste0(save_path, '/pool2_lower600_UMI2000_log_rawCount_Hashtag.csv'))
write.csv(so.p2[['hash.manual']], file = paste0(save_path, '/pool2_lower600_UMI2000_hash_manual.csv'))

### HTO demultiplexing
so.p2 <- HTODemux(so.p2, assay = "hashtag", positive.quantile = 0.99)
write.csv(so.p2[['hashtag_classification']], file = paste0(save_path, '/pool2_lower600_UMI2000_HTOdemux.csv'))

saveRDS(so.p2, file = paste0(save_path, '/pool2_lower600_UMI2000_mt5.rds'))
