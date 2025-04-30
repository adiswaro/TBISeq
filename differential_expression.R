# Differential Expression and Spatial Visualization Functions
# Written by Aditya Swaro, 04/30/2025

getDifferentialExpression <- function(gene) {
  integrated_test <- integrated
  Idents(integrated_test) <- 'conditions'
  
  marker <- FindMarkers(
    integrated_test, ident.1 = 'SHAM', ident.2 = 'TBI',
    features = gene, logfc.threshold = 0, min.pct = 0
  )
  
  marker$Significance <- p.adjust(marker$p_val, method = 'fdr')
  marker <- select(marker, -p_val, -p_val_adj)
  marker$avg_log2FC <- marker$avg_log2FC * (-1)
  
  return(marker)
}

getDifferentialExpressionSpecific <- function(gene, region) {
  integrated_test <- integrated
  integrated_test <- subset(integrated_test, ident = region)
  Idents(integrated_test) <- 'conditions'
  
  marker <- FindMarkers(
    integrated_test, ident.1 = 'SHAM', ident.2 = 'TBI',
    features = gene, logfc.threshold = 0, min.pct = 0,
    recorrect_umi = FALSE
  )
  
  marker$Significance <- p.adjust(marker$p_val, method = 'fdr')
  marker <- select(marker, -p_val, -p_val_adj)
  marker$avg_log2FC <- marker$avg_log2FC * (-1)
  
  return(marker)
}

getSpatialFeature <- function(gene, alpha) {
  shmant <- SpatialFeaturePlot(integrated, features = gene, slot = 'scale.data', images = 'shm1_ant', crop = FALSE, alpha = alpha) +
    ggtitle('SHAM ANTERIOR') + theme_minimal()
  
  shmpos <- SpatialFeaturePlot(integrated, features = gene, slot = 'scale.data', images = 'shm1_pos', crop = FALSE, alpha = alpha) +
    ggtitle('SHAM POSTERIOR') + theme_minimal()
  
  tbiant <- SpatialFeaturePlot(integrated, features = gene, slot = 'scale.data', images = 'tbi2_ant', crop = FALSE, alpha = alpha) +
    ggtitle('TBI ANTERIOR') + theme_minimal()
  
  tbipos <- SpatialFeaturePlot(integrated, features = gene, slot = 'scale.data', images = 'tbi2_pos', crop = FALSE, alpha = alpha) +
    ggtitle('TBI POSTERIOR') + theme_minimal()
  
  return(list(shmant = shmant, shmpos = shmpos, tbiant = tbiant, tbipos = tbipos))
}

getViolinPlot <- function(gene) {
  return(
    VlnPlot(integrated, features = gene, pt.size = 0, split.by = 'conditions', group.by = 'seurat_clusters') +
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "black", colour = NA),
        plot.background = element_rect(fill = "black", colour = "black"),
        panel.border = element_blank(),
        text = element_text(color = "white"),
        axis.text = element_text(color = "white"),
        axis.line = element_line(color = "white"),
        axis.ticks = element_line(color = "white")
      )
  )
}

getSpatialFeatureSpecific <- function(gene, region, alpha) {
  integrated_test <- integrated
  integrated_test$cell.type.stim <- paste(Idents(integrated_test), integrated_test$conditions, sep = "_")
  integrated_test$cell.type <- Idents(integrated_test)
  Idents(integrated_test) <- 'cell.type.stim'
  
  subsetted <- subset(integrated_test, idents = c(paste0(region, '_SHAM'), paste0(region, '_TBI')))
  
  shmant1 <- SpatialFeaturePlot(subsetted, features = gene, pt.size.factor = 0.017 * (1 / integrated@images[["shm1_ant"]]@spot.radius), images = 'shm1_ant', crop = FALSE, alpha = alpha) +
    ggtitle('SHAM ANTERIOR')
  
  shmant2 <- SpatialFeaturePlot(subsetted, features = gene, images = 'shm2_ant', crop = FALSE, alpha = alpha) +
    ggtitle('SHAM ANTERIOR')
  
  shmpos1 <- SpatialFeaturePlot(subsetted, features = gene, images = 'shm1_pos', crop = FALSE, alpha = alpha) +
    ggtitle('SHAM POSTERIOR')
  
  shmpos2 <- SpatialFeaturePlot(subsetted, features = gene, images = 'shm2_pos', crop = FALSE, alpha = alpha) +
    ggtitle('SHAM POSTERIOR')
  
  tbiant1 <- SpatialFeaturePlot(subsetted, features = gene, images = 'tbi1_ant', crop = FALSE, alpha = alpha) +
    ggtitle('TBI ANTERIOR')
  
  tbiant2 <- SpatialFeaturePlot(subsetted, features = gene, images = 'tbi2_ant', crop = FALSE, alpha = alpha) +
    ggtitle('TBI ANTERIOR')
  
  tbipos1 <- SpatialFeaturePlot(subsetted, features = gene, images = 'tbi1_pos', crop = FALSE, alpha = alpha) +
    ggtitle('TBI POSTERIOR')
  
  tbipos2 <- SpatialFeaturePlot(subsetted, features = gene, images = 'tbi2_pos', crop = FALSE, alpha = alpha) +
    ggtitle('TBI POSTERIOR')
  
  return(list(shmant = shmant1, shmpos = shmpos1, tbiant = tbiant1, tbipos = tbipos2))
}

getHeat2 <- function(region) {
  integrated_test <- integrated
  integrated_test$cell.type.stim <- paste(Idents(integrated_test), integrated_test$conditions, sep = "_")
  integrated_test$cell.type <- Idents(integrated_test)
  Idents(integrated_test) <- 'cell.type.stim'
  
  markers <- FindMarkers(integrated_test, ident.1 = paste0(region, '_SHAM'), ident.2 = paste0(region, '_TBI'))
  markers$Significance <- p.adjust(markers$p_val, method = 'fdr')
  markers <- select(markers, -p_val, -p_val_adj)
  markers$Gene <- rownames(markers)
  markers <- select(markers, Gene, everything())
  
  markers <- markers[markers$Significance < 0.05, ]
  markers <- subset(markers, abs(markers$pct.2 - markers$pct.1) > 0.1)
  markers$abs_diff <- abs(markers$pct.2 - markers$pct.1)
  markers <- markers[order(markers$abs_diff), ]
  
  num_to_remove <- floor(0.30 * nrow(markers))
  markers <- markers[(num_to_remove + 1):nrow(markers), ]
  markers$abs_diff <- NULL
  
  markers <- markers[order(markers$Gene, decreasing = TRUE), ]
  markers_long <- reshape2::melt(markers, id.vars = "Gene", measure.vars = c("pct.1", "pct.2"))
  
  return(markers_long)
}
