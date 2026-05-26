source('/home/sjcho/yard/functions/R/FeaturePlot_sjcho.R')
source('/home/sjcho/yard/functions/R/draw_proportion.R')
source('/home/sjcho/yard/functions/R/draw_celltype_wise_dotplot.R')
source('/home/sjcho/yard/functions/R/draw_volcano.R')
source('/home/sjcho/yard/functions/R/draw_stacked_violin.R')
source('/home/sjcho/yard/functions/R/save_ggplot2_to_ppt.R')
source('/home/sjcho/yard/functions/R/draw_violin_split.R')

log_normalize <- function(so, save_path, data_name, pcs, nfeatures = 2000) {
    so <- NormalizeData(object = so)
    so <- FindVariableFeatures(object = so, nfeatures = nfeatures)
    so <- ScaleData(object = so)
    so <- RunPCA(object = so, npcs = 50)
    p <- ElbowPlot(so, ndims = 50, reduction = "pca")
    ggsave(p, file = paste0(save_path, '/', data_name, '_elbowplot.png'), width = 5, height = 5)
    so <- so %>%
        FindNeighbors(reduction = 'pca') %>%
        FindClusters(resolution = 1) 
    so <- RunUMAP(so, dims = 1:pcs, reduction = 'pca')
    return(so)
}

plot_lineage_features <- function(so, save_path, data_name, RdBu = T) {
    p <- fp_sjcho(so, features = c('Epcam', 'Cdh1', 'Ptprc', 'Ptprb', 
                                'Pecam1', 'Pdgfra', 'Col1a1', 'Mki67'), ncol = 4)
    if (RdBu) {
        library(RColorBrewer)
        p <- p & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
    }
    ggsave(p, file = paste0(save_path, data_name, '_lineage_features.png'), width = 13, height = 6)
}

plot_immune_features <- function(so, save_path, data_name, RdBu = T) {
    p <- fp_sjcho(so, features = c('Ptprc', 
                                    'Cd3e', 'Trac', 'Cd4', 'Cd8a', 'Foxp3', 'Trdv4',
                                    'Cd19', 'Cd38',
                                    'Cd68', 'Cd83', 'C1qb'), ncol = 4, order = T)
    if (RdBu) {
        library(RColorBrewer)
        p <- p & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
    }
    ggsave(p, file = paste0(save_path, data_name, '_immune_features.png'), width = 13, height = 9)
}

plot_melanocyte_features <- function(so, save_path, data_name, order = F) {
    p <- fp_sjcho(so, features = c('Mitf', 'Sox10',"Pmel","Dct","Mlana"), ncol = 4, order = order)
    ggsave(p, file = paste0(save_path, data_name, '_melanoma_features.png'), width = 13, height = 6)
}

plot_tech <- function(so, save_path, data_name) {
  so$log10_UMI <- log10(so$nCount_RNA + 1)
  so[['mt.pct']] <- PercentageFeatureSet(so, pattern = "^mt-")
  p <- fp_sjcho(so, features = c('log10_UMI', 'nFeature_RNA', 'mt.pct'), ncol = 2)
  # RdBu
  library(RColorBrewer)
  p <- p & scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
  ggsave(p, file = paste0(save_path, data_name, '_tech.png'), width = 6, height = 5)
}

# hash.cols = c('#bb3838', '#bfaf37', '#67b761', '#3f7fbf', '#bb3fbf')
# names(hash.cols) = c('B16', 'bsAb', 'ICB', 'Combi', 'MC38')

hash.cols = c('#d2d2d2', '#9d9d9d', '#6b6b6b', '#3d3d3d')
names(hash.cols) = c('B16', 'bsAb', 'ICB', 'Combi')

################# colors for annotation
### ann_level1
ann_level2_cols = c('#021da0', '#758eff',
                    'black',
                    '#027161', '#00fcda', '#35961d',
                    '#6200ea', '#af7bf8',
                    '#7d693a',
                    '#ff6f00', '#ece448',
                    '#d50000', '#c99090', '#370202', '#f56f6f')
names(ann_level2_cols) <- c('CD4.naive.like', 'CD4.Th1.like',
                            'Treg.IFN.response',
                            'cTreg', 'fTreg',
                            'Treg.LN.homing', 'cTreg.prolif.', 'fTreg.prolif.',
                            'CD8.naive.like', 
                            'CD8.eff.', 'CD8.eff.prolif.',
                            'CD8.exh.', 'Tpex.Sell-', 'Tpex.Sell+', 'CD8.exh.prolif.')

### ann_level2
ann_level1_cols = c('#0024d6', '#8095ff', 
                    '#000000', '#05dabd', '#6200ea', '#00491e',
                    '#7d693a', '#ff6f00', '#d50000')
names(ann_level1_cols) <- c('CD4.Th1.like', 'CD4.naive.like',
                            'Treg.IFN.response', 'Treg.Ikzf2', 'Treg.Nkg7', 'Treg.Klf2',
                            'CD8.naive.like', 'CD8.effector', 'CD8.exhausted')

############## colors for clone size
# colors_clone_size <- c('gray50', '#28192a', '#0037ff', '#00fffb', '#23a300')
# names(colors_clone_size) <- c('Not Detected', 'Singleton', '2-5', '6-20', '>20')
colors_clone_size <- c('gray50', '#28192a', '#0037ff', '#00fffb', '#23a300', '#ffea00')
names(colors_clone_size) <- c('Not Detected', 'Singleton', '2-5', '6-20', '21-100', '>100')

harmony_on_RNA <- function(so, batch, pc.rna, save_path, sample_name) {
    library(harmony)

    # Run Harmony on RNA
    DefaultAssay(so) <- 'RNA'
    so <- RunHarmony(so, group.by.vars = batch, reduction = 'pca', assay.use = 'RNA', reduction.save = 'harmony.rna')
    so <- RunUMAP(so, dims = 1:pc.rna, reduction = 'harmony.rna', 
            reduction.name = "umap.harmony.rna",
            reduction.key = "harmony.rna.UMAP_")

    so <- so %>%
            FindNeighbors(reduction = "harmony.rna", dims = 1:pc.rna, verbose = T) %>%
            FindClusters(resolution = 1, cluster.name = 'harmony.rna.1', verbose = T)

    p <- DimPlot(so, reduction = 'umap.harmony.rna', group.by = 'harmony.rna.1', label = TRUE)
    ggsave(p, file = paste0(save_path, 'harmnoy_', batch, '_', sample_name, '_umap_rna_cluster.png'), width = 5, height = 4)
    p <- DimPlot(so, reduction = 'umap.harmony.rna', group.by = 'harmony.rna.1', split.by = batch, label = TRUE, ncol = 2)
    ggsave(p, file = paste0(save_path, 'harmnoy_', batch, '_', sample_name, 'split_by_', batch, '_umap_rna_cluster.png'), width = 8, height = 7)

    return(so)
}

set_UMAP <- theme(
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  plot.background = element_rect(fill = "transparent", color = NA),
  panel.background = element_rect(fill = "transparent", color = NA),
  axis.line = element_blank(),
  axis.title = element_blank(),
  plot.title = element_blank())

set_UMAP_featureplot <- theme(
  plot.background = element_rect(fill = "transparent", color = NA),
  panel.background = element_rect(fill = "transparent", color = NA),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  axis.line = element_blank(),
  axis.title = element_blank())

############ add clone size
add_clone_info <- function(so) {
    contig_annotation.p2 <- read.csv('/home/sjcho/projects/4-1BB/DA240060/DA240060_41BB_20250204_pool2/outs/multi/vdj_t/all_contig_annotations.csv')
    contig_annotation.p1 <- read.csv('/home/sjcho/projects/4-1BB/DA240060/DA240060_41BB_20250204_pool1/outs/multi/vdj_t/all_contig_annotations.csv')

    ### merge TCR data into seurat object
    library(scRepertoire)
    so$TCR_clone_size <- 0
    so$CTstrict <- 'Not Yet Assigned'

    ### pool2
    p2.contig <- contig_annotation.p2[paste0('pool2_pool2_', contig_annotation.p2$barcode) %in% colnames(so), ]
    p2.contig$sample = NULL

    combined.TCR <- combineTCR(list(p2.contig), 
                            samples = c('pool2'),
                            removeNA = FALSE, 
                            removeMulti = FALSE, 
                            filterMulti = TRUE)

    so$CTstrict[paste0('pool2_', combined.TCR$pool2$barcode)] <- combined.TCR$pool2$CTstrict

    # set removeMulti to FALSE to keep cell barcode with more than 2 chains
    # set filterMulti to TRUE to isolate the top 2 expressed chains in cell barcdoes with multiple chains

    ## add clone size
    clone_size.p2_strict <- table(combined.TCR$pool2$CTstrict)
    combined.TCR$pool2$clone_size.p2_strict <- as.integer(
    clone_size.p2_strict[combined.TCR$pool2$CTstrict]
    )
    so$TCR_clone_size[paste0('pool2_', combined.TCR$pool2$barcode)] = combined.TCR$pool2$clone_size.p2_strict

    ### pool1
    # Only Remain CD8 T cells which has TCR
    p1.contig <- contig_annotation.p1[paste0('pool1_pool1_', contig_annotation.p1$barcode) %in% colnames(so), ]
    p1.contig$sample = NULL
    combined.TCR <- combineTCR(list(p1.contig), 
                            samples = c('pool1'),
                            removeNA = FALSE, 
                            removeMulti = FALSE, 
                            filterMulti = TRUE)

    so$CTstrict[paste0('pool1_', combined.TCR$pool1$barcode)] <- combined.TCR$pool1$CTstrict

    clone_size.p1_strict <- table(combined.TCR$pool1$CTstrict)
    combined.TCR$pool1$clone_size.p1_strict <- as.integer(
    clone_size.p1_strict[combined.TCR$pool1$CTstrict]
    )
    so$TCR_clone_size[paste0('pool1_', combined.TCR$pool1$barcode)] = combined.TCR$pool1$clone_size.p1_strict
    so$Clone_size_group <- 'Not Detected'
    so$Clone_size_group[so$TCR_clone_size == 0] <- 'Not Detected'
    so$Clone_size_group[so$TCR_clone_size == 1] <- 'Singleton'
    so$Clone_size_group[so$TCR_clone_size >= 2 & so$TCR_clone_size <= 5] <- '2-5'
    so$Clone_size_group[so$TCR_clone_size >= 6 & so$TCR_clone_size <= 20] <- '6-20'
    so$Clone_size_group[so$TCR_clone_size >= 21 & so$TCR_clone_size <= 100] <- '21-100'
    so$Clone_size_group[so$TCR_clone_size > 100] <- '>100'
    so$Clone_size_group <- factor(so$Clone_size_group, levels = c('Not Detected', 'Singleton', '2-5', '6-20', '21-100', '>100'))
    return(so)
}

fp_sjcho <- function(so, features, order = FALSE, pt.size = 0.1, reduction = 'umap', ncol = 1, min.cutoff = NA, max.cutoff = NA, split.by = NULL, cells = NULL, alpha = 1) {
    FeaturePlot(so, features = features, order = order, raster = F, reduction = reduction, ncol = ncol, min.cutoff = min.cutoff, max.cutoff = max.cutoff, split.by = split.by, cells = cells, alpha = alpha,
    cols = c("#e0e0e0", "#b2182b"), pt.size = pt.size) & scale_colour_gradientn(colours = rev(c("#300000", "red","lightgray")))}

set_UMAP <- theme(
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  plot.background = element_rect(fill = "transparent", color = NA),
  panel.background = element_rect(fill = "transparent", color = NA),
  axis.line = element_blank(),
  axis.title = element_blank(),
  plot.title = element_blank())


set_UMAP_featureplot <- theme(
  plot.background = element_rect(fill = "transparent", color = NA),
  panel.background = element_rect(fill = "transparent", color = NA),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  axis.line = element_blank(),
  axis.title = element_blank())

set_plt <- theme(legend.position = "none",
          axis.title.y = element_text(size = rel(1)), 
          axis.text.y = element_text(size = rel(1)), 
          plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"))

set_color_RdBu <- 
  scale_colour_gradientn(
    colours  = rev(brewer.pal(11, "RdBu")),
    limits   = function(x) {
      m <- max(abs(x), na.rm = TRUE)
      c(-m, m)                     
    },
    rescaler = ~ scales::rescale_mid(.x, mid = 0)
  )

set_color_matlab_like2 <- ggplot2::scale_colour_gradientn(
  colours  = rev(colorRamps::matlab.like2(11)), 
  limits   = function(x) {
    m <- max(abs(x), na.rm = TRUE)
    c(-m, m)
  },
  rescaler = ~ scales::rescale_mid(.x, mid = 0)
)

set_color_PiYG <- 
  scale_colour_gradientn(
    colours  = rev(brewer.pal(11, "PiYG")),
    limits   = function(x) {
      m <- max(abs(x), na.rm = TRUE)
      c(-m, m)                     
    },
    rescaler = ~ scales::rescale_mid(.x, mid = 0)
  )

set_color_RdYlBu <- 
  scale_colour_gradientn(
    colours  = rev(brewer.pal(11, "RdYlBu"))
  )

log_normalize <- function(so, save_path, data_name, pcs, nfeatures = 2000) {
  so <- NormalizeData(object = so)
  so <- FindVariableFeatures(object = so, nfeatures = nfeatures)
  so <- ScaleData(object = so)
  so <- RunPCA(object = so, npcs = 50)
  p <- ElbowPlot(so, ndims = 50, reduction = "pca")
  ggsave(p, file = paste0(save_path, '/', data_name, '_elbowplot.png'), width = 5, height = 5)
  so <- so %>%
      FindNeighbors(reduction = 'pca', dims = 1:pcs) %>%
      FindClusters(resolution = 1) 
  so <- RunUMAP(so, dims = 1:pcs, reduction = 'pca')
  return(so)
}

log_normalize_wo_UMAP <- function(so, save_path, data_name, pcs, nfeatures = 2000) {
  so <- NormalizeData(object = so)
  so <- FindVariableFeatures(object = so, nfeatures = nfeatures)
  so <- ScaleData(object = so)
  so <- RunPCA(object = so, npcs = 50)
  p <- ElbowPlot(so, ndims = 50, reduction = "pca")
  ggsave(p, file = paste0(save_path, '/', data_name, '_elbowplot.png'), width = 5, height = 5)
  so <- so %>%
      FindNeighbors(reduction = 'pca', dims = 1:pcs) %>%
      FindClusters(resolution = 1) 
  return(so)
}


log_normalize_patient_specific <- function(so, save_path, data_name, pcs, nfeatures = 2000, patient_key) {
  so[['RNA']] <- split(so[['RNA']], f = so.epi.from.T@meta.data[, patient_key])
  so <- NormalizeData(object = so)
  so <- FindVariableFeatures(object = so, nfeatures = nfeatures)
  so <- ScaleData(object = so)
  so <- RunPCA(object = so, npcs = 50)
  p <- ElbowPlot(so, ndims = 50, reduction = "pca")
  ggsave(p, file = paste0(save_path, '/', data_name, '_elbowplot.png'), width = 5, height = 5)
  so <- so %>%
      FindNeighbors(reduction = 'pca', dims = 1:pcs) %>%
      FindClusters(resolution = 1) 
  so <- RunUMAP(so, dims = 1:pcs, reduction = 'pca')
  return(so)
}

log_normalize_harmony <- function(so, save_path, data_name, pcs, batch_name, nfeatures = 2000) {
  so <- NormalizeData(object = so)
  so <- FindVariableFeatures(object = so, nfeatures = nfeatures)
  so <- ScaleData(object = so)
  so <- RunPCA(object = so, npcs = 50)
  p <- ElbowPlot(so, ndims = 50, reduction = "pca")
  ggsave(p, file = paste0(save_path, '/', data_name, '_elbowplot.png'), width = 5, height = 5)
  so <- RunHarmony(so, batch_name)
  so <- so %>% FindNeighbors(reduction = "harmony", dims = 1:pcs)
  so <- so %>% FindClusters(resolution = 1) 
  so <- RunUMAP(so, dims = 1:pcs, reduction = 'harmony')
  return(so)
}

########### contour
library(patchwork)
########## contour plot from somi Kim from LCB, postech ############

## so : seurat object
## name_of_dimension : name of dimension that contains coordinates
## condition_name : name of metadata that contains condition information / must be factor (require levels!)
## ref_condition : name of reference condition (ex : wt)
## exp_condition : name of experimental condition (ex : ko)
## celltype_name : name of metadata that contains celltype information
## nbin : number of bins for density plot
## linewidth : width of the contour line
## type_color : color for each condition
## celltype_color : color for each celltype
## save_path : path to save the plot
## data_name : name of the data
## width : width of the plot
## height : height of the plot
contour_plot <- function(so, name_of_dimension = 'umap', condition_name, ref_condition, exp_condition, celltype_name, nbin = 10, linewidth = 0.25,
                        type_color, celltype_color, save_path, data_name, width, height) {

  # get umap coordinates, condition, and celltype to build densdf
  x <- so@reductions[name_of_dimension][[1]]@cell.embeddings[, 1]
  y <- so@reductions[name_of_dimension][[1]]@cell.embeddings[, 2]
  condition <- so@meta.data[, condition_name]
  celltype <- so@meta.data[, celltype_name]
  densdf <- data.frame(x = x, y = y, type = condition, celltype = celltype)

  p_ref <- ggplot(densdf, aes(x=x, y=y)) +  
      # First add the point layers
      geom_point(data=densdf[densdf$type != ref_condition,], color = 'lightgray', alpha=1) +
      geom_point(data=densdf[densdf$type == ref_condition,], aes(color = celltype), alpha=0.75) +
      scale_color_manual(values = celltype_color) +
      
      # Then add the contour lines (without fill)
      geom_density_2d(data=densdf[densdf$type == ref_condition,], color = type_color[ref_condition], bins=nbin, linewidth  = linewidth) + 
      
      # Set the plot limits
      xlim(min(densdf$x)-1, max(densdf$x)+1) + 
      ylim(min(densdf$y)-1, max(densdf$y)+1) + 
      theme(legend.position = "none",
          text = element_text(size=15),
          panel.background = element_blank(),
          plot.background = element_blank(),
          panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          strip.background = element_rect(fill="lightgray", colour=NA))
  ggsave(p_ref, file = paste0(save_path,'/',  data_name, '_', ref_condition, '_contour.png'), width = width, height = height)
  ggsave(p_ref, file = paste0(save_path,'/',  data_name, '_', ref_condition, '_contour.pdf'), width = width, height = height)
  ggplot2pptx(p_ref, width, height, paste0(save_path, '/', data_name, '_', ref_condition, '_contour.pptx'))

  p_exp <- ggplot(densdf, aes(x=x, y=y)) +  
      # First add the point layers
      geom_point(data=densdf[densdf$type != exp_condition,], color = 'lightgray', alpha=1) +
      geom_point(data=densdf[densdf$type == exp_condition,], aes(color = celltype), alpha=0.75) +
      scale_color_manual(values = celltype_color) + 
      
      # Then add the contour lines (without fill)
      geom_density_2d(data=densdf[densdf$type == exp_condition,], color = type_color[exp_condition], bins=nbin, linewidth  = linewidth) + 
      
      # Set the plot limits
      xlim(min(densdf$x)-1, max(densdf$x)+1) + 
      ylim(min(densdf$y)-1, max(densdf$y)+1) + 
      theme(legend.position = "none",
          text = element_text(size=15),
          panel.background = element_blank(),
          plot.background = element_blank(),
          panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          strip.background = element_rect(fill="lightgray", colour=NA))
  ggsave(p_exp, file = paste0(save_path, '/', data_name, '_', exp_condition, '_contour.png'), width = width, height = height)
  ggsave(p_exp, file = paste0(save_path, '/', data_name, '_', exp_condition, '_contour.pdf'), width = width, height = height)
  ggplot2pptx(p_exp, width, height, paste0(save_path, '/', data_name, '_', exp_condition, '_contour.pptx'))
}


contour_plot_merge <- function(so,
                               name_of_dimension = 'umap',
                               condition_name,
                               celltype_name,
                               nbin = 10,
                               linewidth = 0.25,
                               celltype_color,
                               condition_color,
                               line_alpha = 1,
                               alpha = 0.75,
                               ncol = NULL,
                               nrow = NULL) 
                               {
  coords <- so@reductions[[name_of_dimension]]@cell.embeddings
  densdf <- data.frame(
    x        = coords[, 1],
    y        = coords[, 2],
    type     = so@meta.data[[condition_name]],
    celltype = so@meta.data[[celltype_name]]
  )
  
  conditions <- levels(so@meta.data[, condition_name])
  plots <- vector("list", length(conditions))
  
  # 전체 데이터 범위 계산 (패딩 추가)
  x_range <- range(densdf$x, na.rm = TRUE)
  y_range <- range(densdf$y, na.rm = TRUE)
  x_padding <- diff(x_range) * 0.1  # 15% 패딩
  y_padding <- diff(y_range) * 0.1  # 15% 패딩
  
  x_limits <- c(x_range[1] - x_padding, x_range[2] + x_padding)
  y_limits <- c(y_range[1] - y_padding, y_range[2] + y_padding)
  
  for (i in seq_along(conditions)) {
    cond <- conditions[i]
    p <- ggplot(densdf, aes(x = x, y = y)) +
      geom_point(data = subset(densdf, type != cond),
                 color = 'lightgray', alpha = alpha, size = 0.5) +
      geom_point(data = subset(densdf, type == cond),
                 aes(color = celltype), alpha = alpha, size = 0.5) +
      scale_color_manual(values = celltype_color) +
      geom_density_2d(data = subset(densdf, type == cond),
                      color = condition_color[[cond]],
                      bins = nbin,
                      linewidth = linewidth,
                      alpha = line_alpha) +
      xlim(x_limits) +
      ylim(y_limits) +
      theme_void(base_size = 15) +
      theme(
        legend.position = "none",
        plot.margin = margin(5, 5, 5, 5),  # 더 큰 마진
        panel.background = element_blank(),
        plot.background = element_blank()
      )
      # + ggtitle(cond)
    
    plots[[i]] <- p
  }
  
  combined <- wrap_plots(plots, ncol = ncol, nrow = nrow)
  return(combined)
}


####### save ggplot2pptx

library(officer)
library(rvg)
library(cowplot)
# save ggplot2 object as ppt
save_ppt <- function(plotObj, fig_save, fig_name) {
    read_pptx() |> # build ppt
    add_slide() |> # add slide
    ph_with( # add image
        dml(ggobj = plotObj), 
        location = ph_location_fullsize() 
    ) |>
    print(paste0(fig_save, fig_name, '.pptx')) # save ppt
}


library(officer)
library(rvg)
ggplot2pptx <- function(p, width, height, filename, fontsize = 7) {
  width_in <- width / 2.54
  height_in <- height / 2.54

  p <- p + theme(
    text             = element_text(size = fontsize),
    axis.text        = element_text(size = fontsize),
    axis.title       = element_text(size = fontsize),
    legend.text      = element_text(size = fontsize),
    legend.title     = element_text(size = fontsize),
    plot.title       = element_text(size = fontsize),
    strip.text       = element_text(size = fontsize),
    plot.background  = element_blank(),
    panel.background = element_blank(),
    legend.background = element_blank(),
    legend.key       = element_blank(),
    strip.background = element_blank()
  )

  editable_graph <- dml(ggobj = p)
  doc <- read_pptx()
  doc <- add_slide(doc)
  doc <- ph_with(x = doc, editable_graph,
                 location = ph_location(width = width, height = height))
  print(doc, target = filename)
}
extract_legend <- function(ggplot_obj) {
  tmp <- ggplot_build(ggplot_obj)
  
  if (!is.null(tmp$plot$scales$get_scales("color"))) {
    color_scale <- tmp$plot$scales$get_scales("color")
    scale_colors <- unique(tmp$data[[1]]$colour)
    scale_labels <- if (!is.null(color_scale$labels)) color_scale$labels else levels(tmp$data[[1]]$group)
    scale_name <- if (!is.null(color_scale$name)) color_scale$name else "color"
  } else {
    scale_colors <- NULL
    scale_labels <- NULL
    scale_name <- NULL
  }
  
  if (!is.null(scale_colors)) {
    dummy_data <- data.frame(
      x = rep(1, length(scale_colors)),
      y = rep(1, length(scale_colors)),
      group = factor(seq_along(scale_colors))
    )
  } else {
    dummy_data <- data.frame(x = 1, y = 1, group = factor(1))
  }
  
  legend_plot <- ggplot() +
    theme_void()
  
  if (!is.null(scale_colors)) {
    legend_plot <- legend_plot +
      geom_point(
        data = dummy_data,
        aes(x = x, y = y, color = group),
        show.legend = TRUE
      ) +
      scale_color_manual(
        values = scale_colors,
        labels = scale_labels,
        name = scale_name
      )
  }
  legend_plot <- legend_plot +
    theme(
      legend.position = "center",
      plot.background = element_blank(),
      panel.background = element_blank(),
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )
    
  return(legend_plot)
}
