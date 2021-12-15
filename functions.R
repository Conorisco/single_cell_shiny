c30 <- c(
   
   "black",
   #"dodgerblue2",#
   "#E31A1C", # red
   "green4", #2
   "#FF7F00", # orange
   "green1",#
   "purple",
   "blue1",#
   "deeppink1",
   "darkorange4",#
   "black",
   "gold1",#
   "darkturquoise",#
   "#6A3D9A", # purpl
   "orchid1",#
   "gray70",
   "maroon",
   "palegreen2",
   "#333333",
   "#CAB2D6", # lt purple
   "#FDBF6F", # lt orange
   "khaki2",#
   #Ω
   "skyblue2",
   "steelblue4",#
   "green1",#
   "yellow4",#
   "yellow3",#
   "#FB9A99", # lt pink
   "brown",
   "#000099",
   "#CC3300"
)

c6 <- c("black",# B cells
        "#E31A1C", # red converting
        "green4", #2 MoMacs
        "#FF7F00", # orange, Monocytes
        "blue1",# ResAMcs
        "grey75" #Unnassigned 
)

#"Monocyte-like SCM"  "DC-like SCM"        "Converting SCM"     "C57BL/6-like LCM"   "BALB/c-like LCM"    "M(IL-4) LCM"        "Naive Intermediate" "Proliferating"      "B Cells" 
c9 <- c(
   
   "#FF7F00", # orange,
   "green4", #2 MoMacs
   
   "#E31A1C", # red converting
   
   "blue1",# 
   "slateblue",
   "darkorchid3",
   "grey60",
   
   "greenyellow",
   "black"# B cells
   
) 

#"Monocyte-like SCM"  "DC-like SCM"        "Converting SCM"     "Naive LCM"          "M(IL-4) LCM"        "Naive Intermediate" "Proliferating"      "B Cells"
c8 <- c(
   
   "#FF7F00", # orange,
   "green4", #2 MoMacs
   "#E31A1C", # red converting
   "blue1",# 
   "darkorchid3",
   "grey60",
   "greenyellow",
   "black"# B cells
   
) 

#"Monocyte-like SCM"  "DC-like SCM"        "Converting SCM"     "LCM"                "Naive Intermediate" "Proliferating"      "B Cells"
c7 <- c(
   
   "#FF7F00", # orange,
   "green4", #2 MoMacs
   "#E31A1C", # red converting
   "blue1",# 
   "grey60",
   "greenyellow",
   "black"# B cells
   
) 

cbin <- c(
   "slateblue", # 1
   "blue1",# 2
   "darkorchid3",#3
   "#FF7F00", # orange, 4
   "green4", # MoMacs
   "#E31A1C"# red converting
) 

c25 <- c(
   "dodgerblue2",#
   "#E31A1C", # red
   "green4", #2
   "#FF7F00", # orange
   "green1",#
   "purple",
   "blue1",#
   "gold1",#
   
   "darkorange4",#
   
   
   "black",
   "deeppink1",
   "#6A3D9A", # purple
   
   "orchid1",#
   
   "gray70",
   "maroon",
   
   "palegreen2",
   "#CAB2D6", # lt purple
   "#FDBF6F", # lt orange
   "khaki2",#
   #Ω
   
   "skyblue2",
   
   "steelblue4",#
   "darkturquoise",#
   "green1",#
   "yellow4",#
   "yellow3",#
   "#FB9A99", # lt pink
   "brown")

c3 <- c("grey77",#
        "#E31A1C", # red
        "blue" #2
)
############################################################################################
cluster_vis <- function(cluster_type, DimRed, SingleCellExperiment_object) {
   stopifnot(!missing(SingleCellExperiment_object))
   stopifnot(!missing(DimRed))
   Sample <- as.factor(SingleCellExperiment_object@colData$Sample)
   
   
   if(cluster_type == "RNA Cluster") {
      Colouring <- as.factor(SingleCellExperiment_object@colData$RNA_Cluster)
      pal <- scale_colour_manual(values = c30)
   } else if (cluster_type == "Cell ID") {
      Colouring <- as.factor(SingleCellExperiment_object@colData$Cell_ID)
      pal <- scale_colour_manual(values = c6)
   } else if (cluster_type == "Regulon Cluster") {
      Colouring <- as.factor(SingleCellExperiment_object@colData$Regulon_Cluster)
      pal <- scale_colour_manual(values = c9)
   } else if (cluster_type == "Regulon Binary Cluster") {
      Colouring <- as.factor(SingleCellExperiment_object@colData$Regulon_binary)
      pal <- scale_colour_manual(values = cbin)
   } else if (cluster_type == "Cluster") {
      Colouring <- as.factor(SingleCellExperiment_object@colData$Cluster)
      pal <- scale_colour_manual(values = c30)
   } else if (cluster_type == "Cell Cycle") {
      Colouring <- as.factor(SingleCellExperiment_object@colData$CellCycle)
      pal <- scale_colour_manual(values = c3)
   } else if (cluster_type == "Sample") {
      Colouring <- as.factor(SingleCellExperiment_object@colData$Sample)
      pal <- scale_colour_manual(values = c30)
   } else if (cluster_type == "Cell Type") {
      Colouring <- as.factor(SingleCellExperiment_object@colData$Cell_Type_LCM_collapse)
      pal <- scale_colour_manual(values = c7)
   } else if (cluster_type == "Cell Type Naive") {
      Colouring <- as.factor(SingleCellExperiment_object@colData$Cell_Type_LCM_naive_collapse)
      pal <- scale_colour_manual(values = c8)
   } else {
      Colouring <- as.factor(SingleCellExperiment_object@colData$label)
      pal <- scale_colour_manual(values = c9)
      print("Error incorrect cluster chosen, options are (in quotes): 'Sample', 'RNA Cluster', 'Regulon Cluster' (default for 'Cluster', or label), 'Cell Type' For Reuglon clusters where all LCM clusters are merged, 'Cell Type Naive' For Regulon clusters where only the 2 naive LCM clusters are merged, 'Cell ID' for manually defined cell types, 'Cell Cycle' for Cell clycle classification)")
   }
   
   # Create a holding tibble to select the correct dimred pass the data on via dpyr later on
   if( DimRed %in%  c(reducedDimNames(SingleCellExperiment_object))){
      holding <- (as_tibble(reducedDim(SingleCellExperiment_object, DimRed)))
   } else {
      holding <- (as_tibble(reducedDim(SingleCellExperiment_object, "UMAP_scenic_auc")))
      print(paste0("Incorrect Reduced Dimention name, using default SCENIC Binary UMAP, acceptable formats: ", paste0(reducedDimNames(sce), collapse = " ")))
   }
   
   # code to draw plot
   holding %>%
      dplyr::mutate("Cluster" = Colouring) %>% 
      dplyr::mutate("Sample" = Sample) %>%
      ggplot(aes(V1, V2, color=Cluster, alpha =0.7)) +
      geom_point(size=0.1,alpha=0.5,aes(colour = Cluster)) +
      #scale_colour_gradientn(colours=terrain.colors(7)) +
      #scale_color_gradientn(colours=c("gray75", "yellow", "red3")) +
      xlab("") + ylab("") + 
      # ggtitle(cluster_type) +
      coord_fixed() +
      theme_classic(base_size=14) +
      guides(colour = guide_legend(title=cluster_type,override.aes = list(size=5, alpha=1))) +
      theme(strip.background = element_blank(),
            #strip.text.x     = element_blank(),
            strip.text.x     = element_text(size=18),
            axis.text.x      = element_blank(),
            axis.text.y      = element_blank(),
            axis.ticks       = element_blank(),
            axis.line        = element_blank(),
            #panel.border     = element_blank(),
            #legend.key.size = unit(1,"line"),
            legend.text=element_text(color=Sample,size=16),
            legend.title = element_text(size = 18),
            plot.title = element_blank(),
            panel.border     = element_rect(colour = "gray50", fill=NA, size=0.25),
            panel.spacing.x=unit(0, "lines"), 
            panel.spacing.y=unit(0,"lines"),
            #legend.position = ("none")) +
            panel.background =  element_blank(), 
            panel.grid.major =  element_blank(),
            panel.grid.minor =  element_blank())  +
      pal +
      facet_wrap(Sample, ncol = 2, nrow = 2)
   #ggsave(path = "gg_save_cluster/", filename = paste(DimRed, cluster_type, ".png",sep="_"), width=8, height=8, dpi=300)
}
########################################################################################
cluster_vis_no_facet <- function(cluster_type, DimRed, SingleCellExperiment_object) {
   stopifnot(!missing(SingleCellExperiment_object))
   stopifnot(!missing(DimRed))
   
   Sample <- as.factor(SingleCellExperiment_object@colData$Sample)
   
   
   if(cluster_type == "RNA Cluster") {
      Colouring <- as.factor(SingleCellExperiment_object@colData$RNA_Cluster)
      pal <- scale_colour_manual(values = c30)
   } else if (cluster_type == "Cell ID") {
      Colouring <- as.factor(SingleCellExperiment_object@colData$Cell_ID)
      pal <- scale_colour_manual(values = c6)
   } else if (cluster_type == "Regulon Cluster") {
      Colouring <- as.factor(SingleCellExperiment_object@colData$Regulon_Cluster)
      pal <- scale_colour_manual(values = c9)
   } else if (cluster_type == "Regulon Binary Cluster") {
      Colouring <- as.factor(SingleCellExperiment_object@colData$Regulon_binary)
      pal <- scale_colour_manual(values = cbin)
   } else if (cluster_type == "Cluster") {
      Colouring <- as.factor(SingleCellExperiment_object@colData$Cluster)
      pal <- scale_colour_manual(values = c30)
   } else if (cluster_type == "Cell Cycle") {
      Colouring <- as.factor(SingleCellExperiment_object@colData$CellCycle)
      pal <- scale_colour_manual(values = c3)
   } else if (cluster_type == "Sample") {
      Colouring <- as.factor(SingleCellExperiment_object@colData$Sample)
      pal <- scale_colour_manual(values = c30)
   } else if (cluster_type == "Cell Type") {
      Colouring <- as.factor(SingleCellExperiment_object@colData$Cell_Type_LCM_collapse)
      pal <- scale_colour_manual(values = c7)
   } else if (cluster_type == "Cell Type Naive") {
      Colouring <- as.factor(SingleCellExperiment_object@colData$Cell_Type_LCM_naive_collapse)
      pal <- scale_colour_manual(values = c8)
   } else if (cluster_type == "Regulon Naive") {
      Colouring <- as.factor(SingleCellExperiment_object@colData$Regulon_Naive)
      pal <- scale_colour_manual(values = c8)
      
   } else {
      Colouring <- as.factor(SingleCellExperiment_object@colData$label)
      pal <- scale_colour_manual(values = c9)
      print("Error incorrect cluster chosen, options are (in quotes): 'Sample', 'RNA Cluster', 'Regulon Cluster' (default for 'Cluster', or label), 'Cell Type' For Reuglon clusters where all LCM clusters are merged, 'Cell Type Naive' For Regulon clusters where only the 2 naive LCM clusters are merged, 'Cell ID' for manually defined cell types, 'Cell Cycle' for Cell clycle classification)")
   }
   
   
   # Create a holding tibble to select the correct dimred pass the data on via dpyr later on
   
   if( DimRed %in%  c(reducedDimNames(SingleCellExperiment_object))){
      holding <- (as_tibble(reducedDim(SingleCellExperiment_object, DimRed)))
   } else {
      holding <- (as_tibble(reducedDim(SingleCellExperiment_object, "UMAP_scenic_auc")))
      print(paste0("Incorrect Reduced Dimention name, using default SCENIC Binary UMAP, acceptable formats: ", paste0(reducedDimNames(sce), collapse = " ")))
   }
   
   
   # code to draw plot
   holding %>%
      dplyr::mutate("Cluster" = Colouring) %>% 
      dplyr::mutate("Sample" = Sample) %>%
      ggplot(aes(V1, V2, color=Cluster)) +
      geom_point(size=0.1,alpha=0.4,aes(colour = Cluster)) +
      #scale_colour_gradientn(colours=terrain.colors(7)) +
      #scale_color_gradientn(colours=c("gray75", "yellow", "red3")) +
      xlab("") + ylab("") + 
      ggtitle(cluster_type) +
      coord_fixed() +
      theme_classic(base_size=14) +
      guides(colour = guide_legend(title=cluster_type,override.aes = list(size=5, alpha=1))) +
      theme(strip.background = element_blank(),
            #strip.text.x     = element_blank(),
            strip.text.x     = element_text(size=18),
            axis.text.x      = element_blank(),
            axis.text.y      = element_blank(),
            axis.ticks       = element_blank(),
            axis.line        = element_blank(),
            #panel.border     = element_blank(),
            #legend.key.size = unit(1,"line"),
            legend.text=element_text(color=Sample,size=16),
            legend.title = element_text(size = 18),
            plot.title = element_blank(),
            #panel.border     = element_rect(colour = "gray50", fill=NA, size=0.25),
            panel.border     = element_blank(),
            panel.spacing.x=unit(0, "lines"), 
            panel.spacing.y=unit(0,"lines"),
            #legend.position = ("none")) +
            panel.background =  element_blank(), 
            panel.grid.major =  element_blank(),
            panel.grid.minor =  element_blank())  +
      pal
   # ggsave(path = "outputs/", filename = paste(DimRed, cluster_type, "no_facet.png",sep="_"), width=8, height=8, dpi=300)
}
#####################################################################################
#####################################################################################
gene_vis <- function(geneName, DimRed, Assay, SingleCellExperiment_object) {
   stopifnot(!missing(SingleCellExperiment_object))
   stopifnot(!missing(DimRed))
   
   #Cluster <- as.factor(SingleCellExperiment_object@colData$Cluster)
   Sample <- as.factor(SingleCellExperiment_object@colData$Sample)
   #Cell_ID <- as.factor(SingleCellExperiment_object@colData$Cell_ID)
   
   if(Assay == "logcounts") {
      GeneExp <- logcounts(SingleCellExperiment_object)[geneName,]
   } else if (Assay == "ALRA") {
      GeneExp <- cpm(SingleCellExperiment_object)[geneName,]
   } else if (Assay == "Regulon_binary") {
      alt_sce <- altExp(sce, 'Regulon_binary')
      GeneExp <- counts(alt_sce)[geneName,]
   } else if (Assay == "Regulon_auc") {
      alt_sce <- altExp(sce, 'Regulon_auc')
      GeneExp <- counts(alt_sce)[geneName,]
   } else if (Assay == 'Regulon_binary_manual') {
      alt_sce <- altExp(sce, 'Regulon_binary_manual')
      GeneExp <- counts(alt_sce)[geneName,]
   } else {
      print("Error incorrect Assay chosen, options are (in quotes): 'logcounts' or 'ALRA' for the imputed data,'Regulon_binary', 'Regulon_auc', or 'Regulon_binary_manual' (for taloired/corrected binerizations, prefered for regulon data) ")
   }
   
   # Create a holding tibble to select the correct dimred pass the data on via dpyr later on
   
   if( DimRed %in%  c(reducedDimNames(SingleCellExperiment_object))){
      holding <- (as_tibble(reducedDim(SingleCellExperiment_object, DimRed)))
   } else {
      holding <- (as_tibble(reducedDim(SingleCellExperiment_object, "UMAP_scenic_auc")))
      print(paste0("Incorrect Reduced Dimention name, using default SCENIC Binary UMAP, acceptable formats: ", paste0(reducedDimNames(sce), collapse = " ")))
   }
   
   
   # code to draw plot
   holding %>%
      dplyr::mutate("Gene Expression" = GeneExp) %>% 
      dplyr::mutate("Sample" = Sample) %>%
      ggplot(aes(V1, V2, color=GeneExp)) +
      geom_point(size=0.25,alpha=0.7,aes(colour = GeneExp)) +
      #scale_colour_gradientn(colours=terrain.colors(7)) +
      scale_color_gradientn(colours=c("gray75", "yellow", "red3")) +
      # scale_color_gradientn(colours=c( "gray80","#FFD998","red1","red3")) +
      xlab("") + ylab("") + 
      ggtitle(geneName) +
      coord_fixed() +
      theme_classic(base_size=14) +
      theme(strip.background = element_blank(),
            #strip.text.x     = element_blank(),
            strip.text.x     = element_text(size=18),
            axis.text.x      = element_blank(),
            axis.text.y      = element_blank(),
            axis.ticks       = element_blank(),
            axis.line        = element_blank(),
            #panel.border     = element_blank(),
            #legend.key.size = unit(1,"line"),
            legend.text=element_text(color=Sample,size=12),
            legend.title = element_text(size = 18),
            plot.title = element_text(size = 18),
            panel.border     = element_rect(colour = "gray50", fill=NA, size=0.25),
            panel.spacing.x=unit(0, "lines"), 
            panel.spacing.y=unit(0,"lines"),
            #legend.position = ("none")) +
            panel.background =  element_blank(), 
            panel.grid.major =  element_blank(),
            panel.grid.minor =  element_blank()) +
      #scale_color_manual(values = c30) +
      
      facet_wrap(Sample, ncol = 2, nrow = 2)
   # ggsave(path = "outputs/", filename = paste(Assay, DimRed, geneName, ".png",sep="_"), width=8, height=8, dpi=300)
}
#####################################################################################
gene_vis_no_facet <- function(geneName, DimRed, Assay, SingleCellExperiment_object) {
   stopifnot(!missing(SingleCellExperiment_object))
   stopifnot(!missing(DimRed))
   
   #Cluster <- as.factor(SingleCellExperiment_object@colData$Cluster)
   Sample <- as.factor(SingleCellExperiment_object@colData$Sample)
   #Cell_ID <- as.factor(SingleCellExperiment_object@colData$Cell_ID)
   
   if(Assay == "logcounts") {
      GeneExp <- logcounts(SingleCellExperiment_object)[geneName,]
   } else if (Assay == "ALRA") {
      GeneExp <- cpm(SingleCellExperiment_object)[geneName,]
   } else if (Assay == "Regulon_binary") {
      alt_sce <- altExp(sce, 'Regulon_binary')
      GeneExp <- counts(alt_sce)[geneName,]
   } else if (Assay == "Regulon_auc") {
      alt_sce <- altExp(sce, 'Regulon_auc')
      GeneExp <- counts(alt_sce)[geneName,]
   } else if (Assay == 'Regulon_binary_manual') {
      alt_sce <- altExp(sce, 'Regulon_binary_manual')
      GeneExp <- counts(alt_sce)[geneName,]
   } else {
      print("Error incorrect Assay chosen, options are (in quotes): 'logcounts' or 'ALRA' for the imputed data,'Regulon_binary', 'Regulon_auc', or 'Regulon_binary_manual' (for taloired/corrected binerizations, prefered for regulon data) ")
   }
   
   # Create a holding tibble to select the correct dimred pass the data on via dpyr later on
   
   if( DimRed %in%  c(reducedDimNames(SingleCellExperiment_object))){
      holding <- (as_tibble(reducedDim(SingleCellExperiment_object, DimRed)))
   } else {
      holding <- (as_tibble(reducedDim(SingleCellExperiment_object, "UMAP_scenic_auc")))
      print(paste0("Incorrect Reduced Dimention name, using default SCENIC Binary UMAP, acceptable formats: ", paste0(reducedDimNames(sce), collapse = " ")))
   }
   
   
   # code to draw plot
   holding %>%
      dplyr::mutate("Gene Expression" = GeneExp) %>% 
      dplyr::mutate("Sample" = Sample) %>%
      ggplot(aes(V1, V2, color=GeneExp)) +
      geom_point(size=0.5,alpha=0.7,aes(colour = GeneExp)) +
      #scale_colour_gradientn(colours=terrain.colors(7)) +
      scale_color_gradientn(colours=c("gray75", "yellow", "red3")) +
      xlab("") + ylab("") + 
      ggtitle(geneName) +
      coord_fixed() +
      theme_classic(base_size=14) +
      theme(strip.background = element_blank(),
            #strip.text.x     = element_blank(),
            strip.text.x     = element_text(size=18),
            axis.text.x      = element_blank(),
            axis.text.y      = element_blank(),
            axis.ticks       = element_blank(),
            axis.line        = element_blank(),
            #panel.border     = element_blank(),
            #legend.key.size = unit(1,"line"),
            legend.text=element_text(color=Sample,size=12),
            legend.title = element_text(size = 18),
            plot.title = element_text(size = 18),
            panel.border     = element_rect(colour = "gray50", fill=NA, size=0.25),
            panel.spacing.x=unit(0, "lines"), 
            panel.spacing.y=unit(0,"lines"),
            #legend.position = ("none")) +
            panel.background =  element_blank(), 
            panel.grid.major =  element_blank(),
            panel.grid.minor =  element_blank()) 
   #scale_color_manual(values = c30) +
   
   # ggsave(path = "outputs/", filename = paste(Assay, DimRed, geneName, "no_facet.png",sep="_"), width=8, height=8, dpi=300)
}

gene_vis_no_facet_small <- function(geneName, DimRed, Assay, SingleCellExperiment_object) {
   stopifnot(!missing(SingleCellExperiment_object))
   stopifnot(!missing(DimRed))
   
   #Cluster <- as.factor(SingleCellExperiment_object@colData$Cluster)
   Sample <- as.factor(SingleCellExperiment_object@colData$Sample)
   #Cell_ID <- as.factor(SingleCellExperiment_object@colData$Cell_ID)
   
   if(Assay == "logcounts") {
      GeneExp <- logcounts(SingleCellExperiment_object)[geneName,]
   } else if (Assay == "ALRA") {
      GeneExp <- cpm(SingleCellExperiment_object)[geneName,]
   } else if (Assay == "Regulon_binary") {
      alt_sce <- altExp(sce, 'Regulon_binary')
      GeneExp <- counts(alt_sce)[geneName,]
   } else if (Assay == "Regulon_auc") {
      alt_sce <- altExp(sce, 'Regulon_auc')
      GeneExp <- counts(alt_sce)[geneName,]
   } else if (Assay == 'Regulon_binary_manual') {
      alt_sce <- altExp(sce, 'Regulon_binary_manual')
      GeneExp <- counts(alt_sce)[geneName,]
   } else {
      print("Error incorrect Assay chosen, options are (in quotes): 'logcounts' or 'ALRA' for the imputed data,'Regulon_binary', 'Regulon_auc', or 'Regulon_binary_manual' (for taloired/corrected binerizations, prefered for regulon data) ")
   }
   
   # Create a holding tibble to select the correct dimred pass the data on via dpyr later on
   
   if( DimRed %in%  c(reducedDimNames(SingleCellExperiment_object))){
      holding <- (as_tibble(reducedDim(SingleCellExperiment_object, DimRed)))
   } else {
      holding <- (as_tibble(reducedDim(SingleCellExperiment_object, "UMAP_scenic_auc")))
      print(paste0("Incorrect Reduced Dimention name, using default SCENIC Binary UMAP, acceptable formats: ", paste0(reducedDimNames(sce), collapse = " ")))
   }
   
   
   # code to draw plot
   holding %>%
      dplyr::mutate("Gene Expression" = GeneExp) %>% 
      dplyr::mutate("Sample" = Sample) %>%
      ggplot(aes(V1, V2, color=GeneExp)) +
      geom_point(size=0.25,alpha=1,aes(colour = GeneExp)) +
      #scale_colour_gradientn(colours=terrain.colors(7)) +
      scale_color_gradientn(colours=c("gray75", "yellow", "red3")) +
      xlab("") + ylab("") + 
      ggtitle(geneName) +
      coord_fixed() +
      theme_classic(base_size=14) +
      theme(strip.background = element_blank(),
            #strip.text.x     = element_blank(),
            strip.text.x     = element_text(size=18),
            axis.text.x      = element_blank(),
            axis.text.y      = element_blank(),
            axis.ticks       = element_blank(),
            axis.line        = element_blank(),
            #panel.border     = element_blank(),
            #legend.key.size = unit(1,"line"),
            # legend.text=element_text(color=Sample,size=12),
            # legend.title = element_text(size = 18),
            plot.title = element_text(size = 20),
            panel.border     = element_rect(colour = "gray25", fill=NA, size=0.6),
            panel.spacing.x=unit(0, "lines"), 
            panel.spacing.y=unit(0,"lines"),
            legend.position = ("none"),
            panel.background =  element_blank(), 
            panel.grid.major =  element_blank(),
            panel.grid.minor =  element_blank(),
            plot.margin = unit(c(0, 0, 0, 0), "cm")) 
   #scale_color_manual(values = c30) +
   
   # ggsave(path = "outputs/", filename = paste(Assay, DimRed, geneName, "no_facet.png",sep="_"), width=8, height=8, dpi=300)
}



##############################################################################################################################################
viz_violin_no_facet <- function(geneName, cluster_type, Assay,  SingleCellExperiment_object) {
   stopifnot(!missing(SingleCellExperiment_object))
   stopifnot(!missing(geneName))
   
   Sample <- as.factor(SingleCellExperiment_object@colData$Sample)
   
   
   if(cluster_type == "RNA Cluster") {
      Colouring <- as.factor(SingleCellExperiment_object@colData$RNA_Cluster)
      pal <- scale_fill_manual(values = c30)
   } else if (cluster_type == "Cell ID") {
      Colouring <- as.factor(SingleCellExperiment_object@colData$Cell_ID)
      pal <- scale_fill_manual(values = c6)
   } else if (cluster_type == "Regulon Cluster") {
      Colouring <- as.factor(SingleCellExperiment_object@colData$Regulon_Cluster)
      pal <- scale_fill_manual(values = c9)
   } else if (cluster_type == "Cluster") {
      Colouring <- as.factor(SingleCellExperiment_object@colData$Cluster)
      pal <- scale_fill_manual(values = c30)
   } else if (cluster_type == "Regulon Binary Cluster") {
      Colouring <- as.factor(SingleCellExperiment_object@colData$Regulon_binary)
      pal <- scale_fill_manual(values = cbin)
   } else if (cluster_type == "Cell Cycle") {
      Colouring <- as.factor(SingleCellExperiment_object@colData$CellCycle)
      pal <- scale_fill_manual(values = c3)
   } else if (cluster_type == "Sample") {
      Colouring <- as.factor(SingleCellExperiment_object@colData$Sample)
      pal <- scale_fill_manual(values = c30)
   } else if (cluster_type == "Cell Type") {
      Colouring <- as.factor(SingleCellExperiment_object@colData$Cell_Type_LCM_collapse)
      pal <- scale_fill_manual(values = c7)
   } else if (cluster_type == "Cell Type Naive") {
      Colouring <- as.factor(SingleCellExperiment_object@colData$Cell_Type_LCM_naive_collapse)
      pal <- scale_fill_manual(values = c8)
   } else {
      Colouring <- as.factor(SingleCellExperiment_object@colData$label)
      pal <- scale_fill_manual(values = c9)
      print("Error incorrect cluster chosen, options are (in quotes): 'Sample', 'RNA Cluster', 'Regulon Cluster' (default for 'Cluster', or label), 'Cell Type' For Reuglon clusters where all LCM clusters are merged, 'Cell Type Naive' For Regulon clusters where only the 2 naive LCM clusters are merged, 'Cell ID' for manually defined cell types, 'Cell Cycle' for Cell clycle classification)")
   }
   
   if(Assay == "logcounts") {
      GeneExp <- logcounts(SingleCellExperiment_object)[geneName,]
      type <-""
      
   } else if (Assay == "ALRA") {
      GeneExp <- cpm(SingleCellExperiment_object)[geneName,]
      type <-""
   } else if (Assay == "Regulon_binary") {
      alt_sce <- altExp(sce, 'Regulon_binary')
      GeneExp <- counts(alt_sce)[geneName,]
      type <-" Regulon"
   } else if (Assay == "Regulon_auc") {
      alt_sce <- altExp(sce, 'Regulon_auc')
      GeneExp <- counts(alt_sce)[geneName,]
      type <-" Regulon"
   } else if (Assay == 'Regulon_binary_manual') {
      alt_sce <- altExp(sce, 'Regulon_binary_manual')
      GeneExp <- counts(alt_sce)[geneName,]
      type <-" Regulon"
   } else {
      GeneExp <- logcounts(SingleCellExperiment_object)[geneName,]
      print("Error incorrect Assay chosen, options are (in quotes): 'logcounts' or 'ALRA' for the imputed data,'Regulon_binary', 'Regulon_auc', or 'Regulon_binary_manual' (for taloired/corrected binerizations, prefered for regulon data) ")
   }
   
   as_tibble(GeneExp) %>%
      dplyr::mutate("Cluster" = Colouring ) %>%
      ggplot(aes(x = Cluster, y = value, fill = Cluster)) +
      ylab("") +
      scale_x_discrete(guide = guide_axis(n.dodge = 2))  +
      scale_y_continuous(breaks = seq(0,1,by=1)) +
      ggtitle(paste0(geneName, type)) +
      geom_violin(trim = T, width = 0.95, scale = "width") +
      theme_classic(base_size=14) +
      
      theme(strip.background = element_blank(),
            axis.text.y  = element_blank(),
            axis.text.x      = element_blank(),
            #panel.border     = element_blank(),
            #legend.key.size = unit(1,"line"),
            legend.title = element_text(size = 22),
            plot.title = element_text(size = 22),
            #legend.position = ("none"),
            panel.background =  element_blank(), 
            panel.grid.major =  element_blank(),
            panel.grid.minor =  element_blank()) +
      pal
   
   # ggsave(path = "outputs/", filename = paste(Assay, DimRed, geneName, "no_facet.png",sep="_"), width=8, height=8, dpi=300)
}

viz_violin<- function(geneName, cluster_type, Assay,  SingleCellExperiment_object) {
   stopifnot(!missing(SingleCellExperiment_object))
   stopifnot(!missing(geneName))
   
   Sample <- as.factor(SingleCellExperiment_object@colData$Sample)
   
   
   if(cluster_type == "RNA Cluster") {
      Colouring <- as.factor(SingleCellExperiment_object@colData$RNA_Cluster)
      pal <- scale_fill_manual(values = c30)
   } else if (cluster_type == "Cell ID") {
      Colouring <- as.factor(SingleCellExperiment_object@colData$Cell_ID)
      pal <- scale_fill_manual(values = c6)
   } else if (cluster_type == "Regulon Cluster") {
      Colouring <- as.factor(SingleCellExperiment_object@colData$Regulon_Cluster)
      pal <- scale_fill_manual(values = c9)
   } else if (cluster_type == "Cluster") {
      Colouring <- as.factor(SingleCellExperiment_object@colData$Cluster)
      pal <- scale_fill_manual(values = c30)
   } else if (cluster_type == "Regulon Binary Cluster") {
      Colouring <- as.factor(SingleCellExperiment_object@colData$Regulon_binary)
      pal <- scale_fill_manual(values = cbin)
   } else if (cluster_type == "Cell Cycle") {
      Colouring <- as.factor(SingleCellExperiment_object@colData$CellCycle)
      pal <- scale_fill_manual(values = c3)
   } else if (cluster_type == "Sample") {
      Colouring <- as.factor(SingleCellExperiment_object@colData$Sample)
      pal <- scale_fill_manual(values = c30)
   } else if (cluster_type == "Cell Type") {
      Colouring <- as.factor(SingleCellExperiment_object@colData$Cell_Type_LCM_collapse)
      pal <- scale_fill_manual(values = c7)
   } else if (cluster_type == "Cell Type Naive") {
      Colouring <- as.factor(SingleCellExperiment_object@colData$Cell_Type_LCM_naive_collapse)
      pal <- scale_fill_manual(values = c8)
   } else {
      Colouring <- as.factor(SingleCellExperiment_object@colData$label)
      pal <- scale_fill_manual(values = c9)
      print("Error incorrect cluster chosen, options are (in quotes): 'Sample', 'RNA Cluster', 'Regulon Cluster' (default for 'Cluster', or label), 'Cell Type' For Reuglon clusters where all LCM clusters are merged, 'Cell Type Naive' For Regulon clusters where only the 2 naive LCM clusters are merged, 'Cell ID' for manually defined cell types, 'Cell Cycle' for Cell clycle classification)")
   }
   
   if(Assay == "logcounts") {
      GeneExp <- logcounts(SingleCellExperiment_object)[geneName,]
      type <-""
      
   } else if (Assay == "ALRA") {
      GeneExp <- cpm(SingleCellExperiment_object)[geneName,]
      type <-""
   } else if (Assay == "Regulon_binary") {
      alt_sce <- altExp(sce, 'Regulon_binary')
      GeneExp <- counts(alt_sce)[geneName,]
      type <-" Regulon"
   } else if (Assay == "Regulon_auc") {
      alt_sce <- altExp(sce, 'Regulon_auc')
      GeneExp <- counts(alt_sce)[geneName,]
      type <-" Regulon"
   } else if (Assay == 'Regulon_binary_manual') {
      alt_sce <- altExp(sce, 'Regulon_binary_manual')
      GeneExp <- counts(alt_sce)[geneName,]
      type <-" Regulon"
   } else {
      GeneExp <- logcounts(SingleCellExperiment_object)[geneName,]
      print("Error incorrect Assay chosen, options are (in quotes): 'logcounts' or 'ALRA' for the imputed data,'Regulon_binary', 'Regulon_auc', or 'Regulon_binary_manual' (for taloired/corrected binerizations, prefered for regulon data) ")
   }
   
   as_tibble(GeneExp) %>%
      dplyr::mutate("Cluster" = Colouring ) %>%
      ggplot(aes(x = Cluster, y = value, fill = Cluster)) +
      ylab("") +
      scale_x_discrete(guide = guide_axis(n.dodge = 2))  +
      scale_y_continuous(breaks = seq(0,1,by=1)) +
      ggtitle(paste0(geneName, type)) +
      geom_violin(trim = T, width = 0.95, scale = "width") +
      theme_classic(base_size=14) +
      
      theme(strip.background = element_blank(),
            axis.text.y  = element_blank(),
            axis.text.x      = element_blank(),
            #panel.border     = element_blank(),
            #legend.key.size = unit(1,"line"),
            legend.title = element_text(size = 22),
            plot.title = element_text(size = 22),
            #legend.position = ("none"),
            panel.background =  element_blank(), 
            panel.grid.major =  element_blank(),
            panel.grid.minor =  element_blank()) +
      pal +
      facet_wrap(Sample, ncol = 4, nrow = 1)
   
   # ggsave(path = "outputs/", filename = paste(Assay, DimRed, geneName, "no_facet.png",sep="_"), width=8, height=8, dpi=300)
}