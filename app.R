library(shiny)
library(shinythemes)
library(shinydashboard)
library(shinyWidgets)

library(Seurat)
library(ggplot2)

# Set data directory and Seurat object path
data.dir <- "/Users/dishasharma/Documents/Coding/Rshiny/CZITest"
rdsObj <- "CZI_Endothelial.rds"

# Read the Seurat object
rds <- readRDS(file.path(data.dir, rdsObj))
rds@meta.data$Site <- factor(rds@meta.data$Site)
head(rds@meta.data$Site) # Checking the Site information
print(colnames(rds@meta.data))

rds$Site_Cluster <- paste(rds$Site,rds$integrated_snn_res.0.075,sep = "_")
rds$Site_CellType <- paste(rds$Site,rds$CellType,sep = "_")

ui <- fluidPage(
  navbarPage("CZI Arterial Cell Atlas", theme = shinytheme("sandstone"),
             tabPanel("UMAP",
                      fluidRow(
                        column(2, img(src = "SM.png", height = 70, width = 100)),
                        column(8, align = "center",h1("CZI Arterial Cell Atlas"))
                      ),
                      fluidRow(
                        column(6,wellPanel(selectInput(inputId = "Site",label = "Choose a Site",choices = unique(rds@meta.data$Site)),
                                           actionButton(inputId = "Site1",label = "Update"))),
                        column(6, wellPanel(textInput(inputId = "Gene", label = "Enter a Gene Name", value = ""),actionButton(inputId = "Gene1", label = "Update"))),
                        br(),
                        
                        column(6, align = "center",h3("Site-wise UMAP and Clusters")),
                        column(6, align = "center",h3("Gene Expression")),
                        
                        column(3,plotOutput("UMAP_SiteWise")),
                        column(3,plotOutput("UMAP_SeuratClusters")),
                        column(3,plotOutput("DotPlot_Gene1")),
                        column(3,plotOutput("FeaturePlot_Gene1")),
                        br(),
                        
                        column(3, downloadButton("download_UMAP_SiteWise", "Download UMAP SiteWise")),
                        column(3, downloadButton("download_UMAP_SeuratClusters", "Download UMAP Seurat Clusters")),
                        column(3, downloadButton("download_DotPlot_Gene1", "Download DotPlot Gene1")),
                        column(3, downloadButton("download_FeaturePlot_Gene1", "Download FeaturePlot Gene1")),
                        br(),
                        
                        column(6,plotOutput("VlnPlot_SiteWiseForGene1")),
                        column(6,plotOutput("VlnPlot_ClusterWiseForGene1")),
                        br(),
                        
                        column(6,downloadButton("download_VlnPlot_SiteWiseForGene1", "Download VlnPlot for Gene-SiteWise")),
                        column(6,downloadButton("download_VlnPlot_ClusterWiseForGene1", "Download VlnPlot for Gene -SeuratClusters")),
                        br()
                      )
             ),
             tabPanel("DEG",
                      fluidRow(
                        column(2, img(src = "SM.png", height = 70, width = 100)),
                        column(8, align = "center",h1("Differential Expression"))
                      ),
                      dashboardPage(
                        dashboardHeader(title = "Differential Expression"),
                        dashboardSidebar(
                          sidebarMenu(
                            menuItem("Site", tabName = "Site", icon = icon("globe")),
                            menuItem("Cluster", tabName = "Cluster", icon = icon("flask")),
                            menuItem("CellType", tabName = "CellType", icon = icon("cogs")),
                            menuItem("Site_Cluster", tabName = "Site_Cluster", icon = icon("users")),
                            menuItem("Site_CellType", tabName = "Site_CellType", icon = icon("database"))
                          )
                        ),
                        dashboardBody(
                          
                          
                          tags$style(HTML("
      /* Set color for menu items (inactive) */
      .main-sidebar .sidebar .sidebar-menu a {
        color: #D2B48C;               /* Text color for inactive menu items */
      }

      /* Set color for active menu item */
      .main-sidebar .sidebar .sidebar-menu .active a {
        background-color: #C9A77D;     /* Background color for active menu item */
        color: white;                  /* Text color for active menu item */
      }

      /* Set color for menu item on hover */
      .main-sidebar .sidebar .sidebar-menu a:hover {
        background-color: #9E7B5E;    /* Hover background color */
        color: white;                 /* Text color on hover */
      }

      /* Icon color in menu items */
      .main-sidebar .sidebar .sidebar-menu .fa {
        color: #A56C42;  /* Set icon color */
      }

      /* Set background color for sidebar when the menu is open */
      .main-sidebar .sidebar {
        background-color: #6F4F27;  /* Sidebar background color */
      }
    ")),
                          tabItems(
                            tabItem(tabName = "Site",
                                    column(6,wellPanel(selectInput("Ident1", "Ident1",choices = sort(unique(rds@meta.data$Site)),
                                                                   multiple = TRUE),actionButton(inputId = "DEGSite",label = "Update"))),
                                    column(6,wellPanel(selectInput("Ident2","Ident2",choices = unique(rds@meta.data$Site),
                                                                   multiple = TRUE),actionButton(inputId = "DEGSite",label = "Update"))),
                                    fluidRow(
                                      column(12, align = "center",h3("Site-Wise Data Table")),
                                      column(12,dataTableOutput("CZI_Site_DEG")),
                                      downloadButton("CZI_Site_DEG_download", "Download Site-Specific Markers Data Table")
                                    )
                            ),
                            tabItem(tabName = "Cluster",
                                    column(6,wellPanel(selectInput("Ident1", "Ident1",choices = sort(unique(rds@meta.data$integrated_snn_res.0.075)),
                                                                   multiple = TRUE))),
                                    column(6,wellPanel(selectInput("Ident2","Ident2",choices = unique(rds@meta.data$integrated_snn_res.0.075),
                                                                   multiple = TRUE),actionButton(inputId = "DEGCluster",label = "Update"))),
                                    fluidRow(
                                      column(12, align = "center",h3("Cluster-Wise Data Table")),
                                      column(12,dataTableOutput("CZI_Cluster_DEG")),
                                      downloadButton("CZI_Cluster_DEG_download", "Download Cluster Markers Data Table")
                                    )),
                            tabItem(tabName = "CellType",
                                    column(6,wellPanel(selectInput("Ident1", "Ident1",choices = sort(unique(rds@meta.data$CellType)),
                                                                   multiple = TRUE))),
                                    column(6,wellPanel(selectInput("Ident2","Ident2",choices = unique(rds@meta.data$CellType),
                                                                   multiple = TRUE),actionButton(inputId = "DEGCellType",label = "Update"))),
                                    fluidRow(
                                      column(12, align = "center",h3("CellType-Wise Data Table")),
                                      column(12,dataTableOutput("CZI_CellType_DEG")),
                                      downloadButton("CZI_CellType_DEG_download", "Download Cell-type Specific Markers Data Table")
                                    )),
                            tabItem(tabName = "Site_Cluster",
                                    column(6,wellPanel(selectInput("Ident1", "Ident1",choices = sort(unique(rds@meta.data$Site_Cluster)),
                                                                   multiple = TRUE))),
                                    column(6,wellPanel(selectInput("Ident2","Ident2",choices = unique(rds@meta.data$Site_Cluster),
                                                                   multiple = TRUE),actionButton(inputId = "DEGSiteCluster",label = "Update"))),
                                    fluidRow(
                                      column(12, align = "center",h3("Site-Cluster-Wise Data Table")),
                                      column(12,dataTableOutput("CZI_SiteCluster_DEG")),
                                      downloadButton("CZI_SiteCluster_DEG_download", "Download Site-Cluster_Specific Markers Data Table")
                                    )),
                            tabItem(tabName = "Site_CellType",
                                    column(6,wellPanel(selectInput("Ident1", "Ident1",choices = sort(unique(rds@meta.data$Site_CellType)),
                                                                   multiple = TRUE))),
                                    column(6,wellPanel(selectInput("Ident2","Ident2",choices = unique(rds@meta.data$Site_CellType),
                                                                   multiple = TRUE),actionButton(inputId = "DEGSiteCellType",label = "Update"))),
                                    fluidRow(
                                      column(12, align = "center",h3("Site-CellType-Wise Data Table")),
                                      column(12,dataTableOutput("CZI_SiteCellType_DEG")),
                                      downloadButton("CZI_SiteCellType_DEG_download", "Download Site-CellType_Specific Markers Data Table")
                                    ))
                          )
                        )
                        
                        
                      )
             ),
             tabPanel("ModuleScore",
                      fluidRow(
                        column(2,img(src = "SM.png", height = 70, width = 70)),
                        column(8, align = "center",h1("Module Score Analysis"))
                      ),
                      fluidRow(
                        column(6,fileInput("upload_GeneScore", "Upload Gene Names"), actionButton(inputId = "upload_GeneScorego",label = "Update")),
                      ),
                      fluidRow(
                        column(6, plotOutput("GeneScore_RNA_VlnPlot")),
                        column(6, plotOutput("GeneScore_RNA_FeaturePlot")),
                      ),
                      fluidRow(
                        column(6,downloadButton("CZI_SiteWise_GeneScore_VlnPlot", "Download VlnPlot for Site-wise Module Score")),
                        column(6,downloadButton("CZI_SiteWise_GeneScore_FeaturePlot", "Download FeaturePlot for Cells-wise Module Score")),
                      )
             )
             
  )
  
)
server <- function(input, output) {
  rds_Site <- reactive({
    req(input$Site)  # Ensure input$Site is available before proceeding
    print(input$Site) # Convert 'Site' to character for comparison if necessary
    rds@meta.data$Site <- as.character(rds@meta.data$Site)
    subset(rds, subset = Site %in% input$Site)  # Subset based on the selected Site
  })
  
  # Site-wise UMAP plot
  Site_wiseUMAP <- eventReactive(input$Site1, {
    seurat_obj <- rds_Site()  # Get the Seurat object based on the selected site
    DefaultAssay(seurat_obj) <- "RNA"
    Idents(seurat_obj) <- 'Site'
    DimPlot(seurat_obj, reduction = "umap") + ggtitle(paste(input$Site))
  })
  output$UMAP_SiteWise <- renderPlot({
    Site_wiseUMAP()
  })
  output$download_UMAP_SiteWise <- downloadHandler(
    filename = function() { paste(input$Site, "UMAP.pdf", sep = '_') },
    content = function(file) {
      ggsave(file, plot = Site_wiseUMAP(), device = "pdf", width = 5, height = 5)
    }
  )
  
  # Site-wise Seurat Clusters UMAP plot
  Site_SeuratClustersUMAP <- eventReactive(input$Site1, {
    seurat_obj <- rds_Site()  # Get the Seurat object based on the selected site
    DefaultAssay(seurat_obj) <- "RNA"
    Idents(seurat_obj) <- 'integrated_snn_res.0.075'
    DimPlot(seurat_obj, reduction = "umap", group.by = "integrated_snn_res.0.075") + ggtitle(paste("Seurat Clusters"))
  })
  output$UMAP_SeuratClusters <- renderPlot({
    Site_SeuratClustersUMAP()
  })
  output$download_UMAP_SeuratClusters <- downloadHandler(
    filename = function() { paste(input$Site, "Seurat_clusters_UMAP.pdf", sep = '_') },
    content = function(file) {
      ggsave(file, plot = Site_SeuratClustersUMAP(), device = "pdf", width = 5, height = 5)
    }
  )
  
  # Gene-wise DotPlot
  Gene_wiseDotPlot <- eventReactive(input$Gene1, {
    req(input$Gene)  # Ensure input$Gene is not empty
    if (input$Gene != "") {
      seurat_obj <- rds  # Get the Seurat object based on the selected site
      DefaultAssay(seurat_obj) <- "RNA"
      Idents(seurat_obj) <- 'Site'
      DotPlot(object = seurat_obj, features = input$Gene) + ggtitle(paste("Gene Expression for:", input$Gene))
    } else {
      return(NULL)
    }
  })
  output$DotPlot_Gene1 <- renderPlot({
    Gene_wiseDotPlot()
  })
  output$download_DotPlot_Gene1 <- downloadHandler(
    filename = function() { paste(input$Gene, "Sitewise_DotPlot.pdf", sep = '_') },
    content = function(file) {
      ggsave(file, plot = Gene_wiseDotPlot(), device = "pdf", width = 5, height = 5)
    }
  )
  
  
  # Gene-wise FeaturePlot
  Gene_wiseFeaturePlot <- eventReactive(input$Gene1, {
    req(input$Gene)  # Ensure input$Gene is not empty
    if (input$Gene != "") {
      seurat_obj <- rds_Site()  # Get the Seurat object based on the selected site
      DefaultAssay(seurat_obj) <- "RNA"
      Idents(seurat_obj) <- 'Site'
      FeaturePlot(object = seurat_obj, features = input$Gene) + ggtitle(paste("Feature Plot for Gene:", input$Gene))
    } else {
      return(NULL)
    }
  })
  output$FeaturePlot_Gene1 <- renderPlot({
    Gene_wiseFeaturePlot()
  })
  output$download_FeaturePlot_Gene1 <- downloadHandler(
    filename = function() { paste(input$Gene, "Sitewise_FeaturePlot.pdf", sep = '_') },
    content = function(file) {
      ggsave(file, plot = Gene_wiseFeaturePlot(), device = "pdf", width = 5, height = 5)
    }
  )
  
  # VlnPlot for Gene-SiteWise
  Gene_SitewiseVlnPlot <- eventReactive(input$Gene1, {
    req(input$Gene)  # Ensure input$Gene is not empty
    if (input$Gene != "") {
      seurat_obj <- rds  # Get the Seurat object based on the selected site
      DefaultAssay(seurat_obj) <- "RNA"
      Idents(seurat_obj) <- 'Site'
      VlnPlot(object = seurat_obj, features = input$Gene) + ggtitle(paste("Violin Plot for Gene:", input$Gene))
    } else {
      return(NULL)
    }
  })
  output$VlnPlot_SiteWiseForGene1 <- renderPlot({
    Gene_SitewiseVlnPlot()
  })
  output$download_VlnPlot_SiteWiseForGene1 <- downloadHandler(
    filename = function() { paste(input$Gene, "Sitewise_VlnPlot.pdf", sep = '_') },
    content = function(file) {
      ggsave(file, plot = Gene_SitewiseVlnPlot(), device = "pdf", width = 5, height = 5)
    }
  )
  
  # VlnPlot for Gene-ClusterWise
  Gene_ClusterwiseVlnPlot <- eventReactive(input$Gene1, {
    req(input$Gene)  # Ensure input$Gene is not empty
    if (input$Gene != "") {
      seurat_obj <- rds_Site()  # Get the Seurat object based on the selected site
      DefaultAssay(seurat_obj) <- "RNA"
      Idents(seurat_obj) <- 'integrated_snn_res.0.075'
      VlnPlot(object = seurat_obj, features = input$Gene) + ggtitle(paste("Violin Plot for Gene:", input$Gene))
    } else {
      return(NULL)
    }
  })
  output$VlnPlot_ClusterWiseForGene1 <- renderPlot({
    Gene_ClusterwiseVlnPlot()
  })
  output$download_VlnPlot_ClusterWiseForGene1 <- downloadHandler(
    filename = function() { paste(input$Gene, "Clusterwise_VlnPlot.pdf", sep = '_') },
    content = function(file) {
      ggsave(file, plot = Gene_ClusterwiseVlnPlot(), device = "pdf", width = 5, height = 5)
    }
  )
  
  ########## DEG ##############
  rds_DEG_Site <- eventReactive(input$DEGSite, {
    DefaultAssay(rds) <- "RNA"
    Site_FM <- FindMarkers(rds, ident.1 = input$Ident1, ident.2 = input$Ident2, group.by = "Site", min.pct = 0, logfc.threshold = 0.1)
    Site_FM$Gene <- rownames(Site_FM)
    Site_FM
  })
  
  output$CZI_Site_DEG <- renderDataTable({
    rds_DEG_Site()
  })
  
  output$CZI_Site_DEG_download <- downloadHandler(
    filename = function(){paste(input$Ident1,"vs",input$Ident2,"SiteWise_markers.csv", sep='_')},
    content = function(fname){
      write.csv(rds_DEG_Site(), fname)
    }
  )
  
  
  rds_Cluster <- eventReactive(input$DEGCluster, {
    DefaultAssay(rds) <- "RNA"
    Cluster_FM <- FindMarkers(rds, ident.1 = input$Ident1, ident.2 = input$Ident2, group.by = "integrated_snn_res.0.075", min.pct = 0, logfc.threshold = 0.1)
    Cluster_FM$Gene <- rownames(Cluster_FM)
    Cluster_FM
  })
  
  output$CZI_Cluster_DEG <- renderDataTable({
    rds_Cluster()
  })
  
  output$CZI_Cluster_DEG_download <- downloadHandler(
    filename = function(){paste(input$Ident1,"vs",input$Ident2,"ClusterWise_markers.csv", sep='_')},
    content = function(fname){
      write.csv(rds_Cluster(), fname)
    }
  )
  
  
  rds_CellType <- eventReactive(input$DEGCellType, {
    DefaultAssay(rds) <- "RNA"
    CellType_FM <- FindMarkers(rds, ident.1 = input$Ident1, ident.2 = input$Ident2, group.by = "CellType", min.pct = 0, logfc.threshold = 0.1)
    CellType_FM$Gene <- rownames(CellType_FM)
    CellType_FM
  })
  
  output$CZI_CellType_DEG <- renderDataTable({
    rds_CellType()
  })
  
  output$CZI_CellType_DEG_download <- downloadHandler(
    filename = function(){paste(input$Ident1,"vs",input$Ident2,"CellType_markers.csv", sep='_')},
    content = function(fname){
      write.csv(rds_CellType(), fname)
    }
  )
  
  
  rds_Site_Cluster <- eventReactive(input$DEGSiteCluster, {
    DefaultAssay(rds) <- "RNA"
    Site_Cluster_FM <- FindMarkers(rds, ident.1 = input$Ident1, ident.2 = input$Ident2, group.by = "Site_Cluster", min.pct = 0, logfc.threshold = 0.1)
    Site_Cluster_FM$Gene <- rownames(Site_Cluster_FM)
    Site_Cluster_FM
  })
  
  output$CZI_SiteCluster_DEG <- renderDataTable({
    rds_Site_Cluster()
  })
  
  output$CZI_SiteCluster_DEG_download <- downloadHandler(
    filename = function(){paste(input$Ident1,"vs",input$Ident2,"Site_Cluster_markers.csv", sep='_')},
    content = function(fname){
      write.csv(rds_Site_Cluster(), fname)
    }
  )
  
  
  rds_Site_CellType <- eventReactive(input$DEGSiteCellType, {
    DefaultAssay(rds) <- "RNA"
    Site_CellType_FM <- FindMarkers(rds, ident.1 = input$Ident1, ident.2 = input$Ident2, group.by = "Site_CellType", min.pct = 0, logfc.threshold = 0.1)
    Site_CellType_FM$Gene <- rownames(Site_CellType_FM)
    Site_CellType_FM
  })
  
  output$CZI_SiteCellType_DEG <- renderDataTable({
    rds_Site_CellType()
  })
  
  output$CZI_SiteCellType_DEG_download <- downloadHandler(
    filename = function(){paste(input$Ident1,"vs",input$Ident2,"Site_Site_CellType_markers.csv", sep='_')},
    content = function(fname){
      write.csv(rds_Site_CellType(), fname)
    }
  )
  
  ### GeneScore Analysis ####
  # Render the violin plot for the uploaded gene scores
  GeneScore_analysis <- eventReactive(input$upload_GeneScore, {
    signature.names <- read.csv(input$upload_GeneScore$datapath, header = TRUE)  # Adjust for file format
    signature.names <- as.character(signature.names$Gene)  # Assuming the file has a "Gene" column
    
    DefaultAssay(rds) <- "RNA"
    rds <- AddModuleScore(rds, features = signature.names, ctrl = 5, name = 'ModuleScore1')
    VlnPlot(rds, features = "ModuleScore1", group.by = "Site") + 
      ggtitle("Module Score by Site")
  })
  
  output$GeneScore_RNA_VlnPlot <- renderPlot({
    GeneScore_analysis()
  })
  
  output$CZI_SiteWise_GeneScore_VlnPlot <- downloadHandler(
    filename = function() { paste("ModuleScore_VlnPlot.png", sep='_') },
    content = function(fname) {
      ggsave(fname, plot = GeneScore_analysis(), device = "png", width = 5, height = 5)
    }
  )
  
  GeneScore_analysis_FP <- eventReactive(input$upload_GeneScore, {
    signature.names <- read.csv(input$upload_GeneScore$datapath, header = TRUE) 
    signature.names <- as.character(signature.names$Gene)
    DefaultAssay(rds) <- "RNA"
    rds <- AddModuleScore(rds, features = signature.names, ctrl = 5, name = 'ModuleScore1')
    FeaturePlot(rds, features = "ModuleScore1") + 
      ggtitle("Module Score across Cells")
  })
  output$GeneScore_RNA_FeaturePlot <- renderPlot({
    GeneScore_analysis_FP()
  })
  output$CZI_SiteWise_GeneScore_FeaturePlot <- downloadHandler(
    filename = function() { paste("ModuleScore_FeaturePlot.png", sep='_') },
    content = function(fname) {
      ggsave(fname, plot = GeneScore_analysis_FP(), device = "png", width = 5, height = 5)
    }
  )
  
  
}

shinyApp(ui = ui, server = server)