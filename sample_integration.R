#Integration of all samples. This script requires all 8 samples to be preloaded
#Written by Angela Zhang and Aditya Swaro, 04/30/2025


library(tidyverse)

#After loading slices---------------------------------------------------------

#Setting orig.ident (labels all imagse correctly)-----------------------------

S1$orig.ident <- 'shm1_ant'
S2$orig.ident <- 'shm1_pos'
S3$orig.ident <- 'tbi1_ant'
S4$orig.ident <- 'tbi1_pos'
S5$orig.ident <- 'shm2_ant'
S6$orig.ident <- 'shm2_pos'
S7$orig.ident <- 'tbi2_ant'
S8$orig.ident <- 'tbi2_pos'

#Quality Control - Removing MT Genes -----------------------------------------

S1 <- S1[!grepl("^mt-", rownames(S1)), ]
S2 <- S2[!grepl("^mt-", rownames(S2)), ]
S3 <- S3[!grepl("^mt-", rownames(S3)), ]
S4 <- S4[!grepl("^mt-", rownames(S4)), ]
S5 <- S5[!grepl("^mt-", rownames(S5)), ]
S6 <- S6[!grepl("^mt-", rownames(S6)), ]
S7 <- S7[!grepl("^mt-", rownames(S7)), ]
S8 <- S8[!grepl("^mt-", rownames(S8)), ]

#SCTransform Normalization ---------------------------------------------------

S1[["RNA"]] <- S1[["Spatial"]]
S2[["RNA"]] <- S2[["Spatial"]]
S3[["RNA"]] <- S3[["Spatial"]]
S4[["RNA"]] <- S4[["Spatial"]]
S5[["RNA"]] <- S5[["Spatial"]]
S6[["RNA"]] <- S6[["Spatial"]]
S7[["RNA"]] <- S7[["Spatial"]]
S8[["RNA"]] <- S8[["Spatial"]]


S1 <- SCTransform(S1, assay = 'RNA')
S2 <- SCTransform(S2, assay = 'RNA')
S3 <- SCTransform(S3, assay = 'RNA')
S4 <- SCTransform(S4, assay = 'RNA')
S5 <- SCTransform(S5, assay = 'RNA')
S6 <- SCTransform(S6, assay = 'RNA')
S7 <- SCTransform(S7, assay = 'RNA')
S8 <- SCTransform(S8, assay = 'RNA')


#Integrating slices ----------------------------------------------------------

IntegrateVisiumSlices <- function(list.slices) {
  options(future.globals.maxSize = 20000 * 1024^2)
  st.list <- list(S1, S2, S3, S4, S5, S6, S7, S8)
  # st.list <- lapply(st.list, SCTransform,
  #                   assay = "RNA", method = "poisson")
  st.features <- SelectIntegrationFeatures(st.list, nfeatures = 3000,
                                           verbose = FALSE)
  st.list <- PrepSCTIntegration(object.list = st.list,
                                anchor.features = st.features,
                                verbose = T)
  int.anchors <- FindIntegrationAnchors(object.list = st.list,
                                        normalization.method = "SCT",
                                        anchor.features = st.features)
  integrated <- IntegrateData(anchorset = int.anchors,
                              normalization.method = "SCT",
                              verbose = T)
  
}

integrated_slices <- IntegrateVisiumSlices(list(S1, S2, S3, S4, S5, S6, S7, S8))

#Performing Dimensional Reduction ---------------------------------------------

RunDimensionalReduction <- function(integrated){
  integrated2 <- RunPCA(integrated)
  integrated2 <- FindNeighbors(integrated2, dims = 1:30)
  integrated2 <- FindClusters(integrated2, resolution = 1.3)
  integrated2 <- RunUMAP(integrated2, dims = 1:30)
}

integrated_slices2 <- integrated

integrated_slices <- RunDimensionalReduction(integrated_slices2)

####CREATING A COLUMN THAT SPECIFIES SHAM OR TBI

integrated_slices@meta.data <- integrated_slices@meta.data |>
  mutate(conditions = case_when(
    orig.ident %in% c("shm1_ant", "shm2_ant", "shm1_pos", "shm2_pos") ~ "SHAM",
    orig.ident %in% c("tbi1_ant", "tbi2_ant", "tbi1_pos", "tbi2_pos") ~ "TBI",
    TRUE ~ NA_character_))


#ROTATING IMAGES

rotimat=function(foo,rotation){
  if(!is.matrix(foo)){
    cat("Input is not a matrix")
    return(foo)
  }
  if(!(rotation %in% c("180","Hf","Vf", "R90", "L90"))){
    cat("Rotation should be either L90, R90, 180, Hf or Vf\n")
    return(foo)
  }
  if(rotation == "180"){
    foo <- foo %>% 
      .[, dim(.)[2]:1] %>%
      .[dim(.)[1]:1, ]
  }
  if(rotation == "Hf"){
    foo <- foo %>%
      .[, dim(.)[2]:1]
  }
  
  if(rotation == "Vf"){
    foo <- foo %>%
      .[dim(.)[1]:1, ]
  }
  if(rotation == "L90"){
    foo = t(foo)
    foo <- foo %>%
      .[dim(.)[1]:1, ]
  }
  if(rotation == "R90"){
    foo = t(foo)
    foo <- foo %>%
      .[, dim(.)[2]:1]
  }
  return(foo)
}

rotateSeuratImage = function(seuratVisumObject, slide = "slice1", rotation="Vf"){
  if(!(rotation %in% c("180","Hf","Vf", "L90", "R90"))){
    cat("Rotation should be either 180, L90, R90, Hf or Vf\n")
    return(NULL)
  }else{
    seurat.visium = seuratVisumObject
    ori.array = (seurat.visium@images)[[slide]]@image
    img.dim = dim(ori.array)[1:2]/(seurat.visium@images)[[slide]]@scale.factors$lowres
    new.mx <- c()  
    # transform the image array
    for (rgb_idx in 1:3){
      each.mx <- ori.array[,,rgb_idx]
      each.mx.trans <- rotimat(each.mx, rotation)
      new.mx <- c(new.mx, list(each.mx.trans))
    }
    
    # construct new rgb image array
    new.X.dim <- dim(each.mx.trans)[1]
    new.Y.dim <- dim(each.mx.trans)[2]
    new.array <- array(c(new.mx[[1]],
                         new.mx[[2]],
                         new.mx[[3]]), 
                       dim = c(new.X.dim, new.Y.dim, 3))
    
    #swap old image with new image
    seurat.visium@images[[slide]]@image <- new.array
    
    ## step4: change the tissue pixel-spot index
    img.index <- (seurat.visium@images)[[slide]]@coordinates
    
    #swap index
    if(rotation == "Hf"){
      seurat.visium@images[[slide]]@coordinates$imagecol <- img.dim[2]-img.index$imagecol
    }
    
    if(rotation == "Vf"){
      seurat.visium@images[[slide]]@coordinates$imagerow <- img.dim[1]-img.index$imagerow
    }
    
    if(rotation == "180"){
      seurat.visium@images[[slide]]@coordinates$imagerow <- img.dim[1]-img.index$imagerow
      seurat.visium@images[[slide]]@coordinates$imagecol <- img.dim[2]-img.index$imagecol
    }
    
    if(rotation == "L90"){
      seurat.visium@images[[slide]]@coordinates$imagerow <- img.dim[2]-img.index$imagecol
      seurat.visium@images[[slide]]@coordinates$imagecol <- img.index$imagerow
    }
    
    if(rotation == "R90"){
      seurat.visium@images[[slide]]@coordinates$imagerow <- img.index$imagecol
      seurat.visium@images[[slide]]@coordinates$imagecol <- img.dim[1]-img.index$imagerow
    }
    
    return(seurat.visium)
  }  
}


integrated_slices2 <- rotateSeuratImage(integrated_slices2,rotation = "Vf", slide = "shm1_ant")
integrated_slices2 <- rotateSeuratImage(integrated_slices2,rotation = "Vf", slide = "shm1_pos")
integrated_slices2 <- rotateSeuratImage(integrated_slices2,rotation = "Vf", slide = "shm2_ant")
integrated_slices2 <- rotateSeuratImage(integrated_slices2,rotation = "Vf", slide = "tbi2_ant")
integrated_slices2 <- rotateSeuratImage(integrated_slices2,rotation = "Vf", slide = "tbi2_pos")
integrated_slices2 <- rotateSeuratImage(integrated_slices2,rotation = "Vf", slide = "tbi1_ant")
integrated_slices2 <- rotateSeuratImage(integrated_slices2,rotation = "Vf", slide = "tbi1_pos")

####CREATING A COLUMN THAT SPECIFIES ANT OR POS


integrated_slices@meta.data <- integrated_slices@meta.data |>
  mutate(slice_conditions = case_when(
    orig.ident %in% c("shm1_ant", "shm2_ant") ~ "Anterior",
    orig.ident %in% c("shm1_pos", "shm2_pos") ~ "Posterior",
    orig.ident %in% c("tbi1_ant", "tbi2_ant") ~ "Anterior",
    orig.ident %in% c("tbi1_pos", "tbi2_pos") ~ "Posterior",
    TRUE ~ NA_character_))

###CREATING A COLUMN WITH READABLE CLUSTER NAMES

new_names <- c("Area 1", "Area 2", "Area 3", "Area 4",
               "Area 5", "Area 6", "Area 7", "Area 8", "Area 9",
               "Area 10", "Area 11", "Area 12", "Area 13", "Area 14",
               "Area 15", "Area 16", "Area 17", "Area 18", "Area 19",
               "Area 20", "Area 21", "Area 22", "Area 23", "Area 24",
               "Area 25", "Area 26", "Area 27", "Area 28", "Area 29",
               "Area 30", "Area 31", "Area 32", "Area 33", "Area 34",
               "Area 35", "Area 36")
names(new_names) <- levels(integrated_slices)

integrated_slices <- RenameIdents(object = integrated_slices, new_names)

integrated_slices@meta.data$seurat_clusters <- recode(integrated_slices@meta.data$seurat_clusters,
                                                      "0" = "Area 1",
                                                      "1" = "Area 2",
                                                      "2" = "Area 3",
                                                      "3" = "Area 4",
                                                      "4" = "Area 5",
                                                      "5" = "Area 6",
                                                      "6" = "Area 7",
                                                      "7" = "Area 8",
                                                      "8" = "Area 9",
                                                      "9" = "Area 10",
                                                      "10" = "Area 11",
                                                      "11" = "Area 12",
                                                      "12" = "Area 13",
                                                      "13" = "Area 14",
                                                      "14" = "Area 15",
                                                      "15" = "Area 16",
                                                      "16" = "Area 17",
                                                      "17" = "Area 18",
                                                      "18" = "Area 19",
                                                      "19" = "Area 20",
                                                      "20" = "Area 21",
                                                      "21" = "Area 22",
                                                      "22" = "Area 23",
                                                      "23" = "Area 24",
                                                      "24" = "Area 25",
                                                      "25" = "Area 26",
                                                      "26" = "Area 27",
                                                      "27" = "Area 28",
                                                      "28" = "Area 29",
                                                      "29" = "Area 30",
                                                      "30" = "Area 31",
                                                      "31" = "Area 32",
                                                      "32" = "Area 33",
                                                      "33" = "Area 34",
                                                      "34" = "Area 35",
                                                      "35" = "Area 36",
)

###SETTING DEFAULT ASSAY TO SCT FOR ALL DOWNSTREAM ANALYSES ON INTEGRATED OBJECT

DefaultAssay(integrated_slices) <- "SCT"

### NEED TO DO THIS STEP BELOW TO RUN FINDMARKERS ON AN SCT ASSAY

integrated_slices <- PrepSCTFindMarkers(integrated_slices, assay = "SCT", verbose = TRUE)

#Exporting integrated object

save(integrated_slices, file = "integrated_v1.RData")