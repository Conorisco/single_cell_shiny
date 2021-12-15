
#Load the required packages
library(shiny)
library(tidyverse)
library(scater)

Sys.setenv(R_MAX_VSIZE = 16e9)
#define colours

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

#Load the data

load(file = "sce.Rdata")

# Data for dropdown list

all_geneNames <- rownames(sce)

all_regulonNames <- rownames(altExp(sce, 'Regulon_auc'))

all_DimReds <- c("UMAP_scenic_auc", "PCA")

all_Assays <- c('logcounts',"Regulon_auc")


# Define UI for application 
ui <- fluidPage(
    # Add some HTML text to the app
    h1("Allen lab Visualisation App"),
    p("Shiny App developed by Conor Finlay: https://github.com/Conorisco"),
    
    
    # Application title
    titlePanel("Gene Expression Plots"),
    
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
          selectInput("regulonName",
                      "Regulon",
                      all_regulonNames),
      #      selectInput("DimRed_input", 
  #                      "Dimention Reduction",
   #                     all_DimReds),
            submitButton("Generate Plot")
        ),
        
        # Show a plot of the generated distribution
        mainPanel(
            plotOutput(outputId = "Sample_plot")
        )
    ),
    #Add a download button
    downloadButton('foo')
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    
    output$Sample_plot<- renderPlot({
        # use input$x from ui.R to provide input to the plot
       
        #geneName <- "Chil3"
      alt_sce <- altExp(sce, 'Regulon_auc')
      GeneExp <- counts(alt_sce)[input$regulonName,]
        
        Sample <- as.factor(sce@colData$Sample)
        holding <- (as_tibble(reducedDim(sce, "UMAP_scenic_auc")))
        
           
      #  if(input$DimRed %in%  c(reducedDimNames(sce))){
      #      holding <- (as_tibble(reducedDim(sce, input$DimRed)))
      #  } else {
      #      holding <- (as_tibble(reducedDim(sce, "UMAP_scenic_auc")))
      #      print(paste0("Incorrect Reduced Dimention name, using default SCENIC Binary UMAP, acceptable formats: ", paste0(reducedDimNames(sce), collapse = " ")))
     #   }
      
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
            ggtitle(input$geneName) +
            coord_fixed() +
            theme_classic(base_size=24) +
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
    }, 
    
    
    height = 800)
}

# Run the application 
shinyApp(ui = ui, server = server)
