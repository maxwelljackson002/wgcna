library(shiny)
library(openxlsx)
library(WGCNA)
library(igraph)
library(BiocManager)
library(tidyverse)
options(repos = BiocManager::repositories())

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("Weighted Gene Co-Expression Network Analysis"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      fileInput("dataFile", "Upload Gene Expression Data", accept = c(".csv", ".xlsx")),
#      sliderInput("bins",
 #                 "Number of bins:",
  #                min = 1,
   #               max = 2700,
    #              value = 30),
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(
        tabPanel("Module Assignment", dataTableOutput("moduleTable")),
        tabPanel("Module Dendrogram", plotOutput("moduleDendrogram")),
        tabPanel("Module Network", plotOutput("moduleNetwork")),
#       tabPanel("Module Detection", plotOutput("moduleDetectionTable"))
      )
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  # Process the uploaded file and display as a table
  output$moduleTable <- renderDataTable({
    # Check if a file is uploaded
    if (!is.null(input$dataFile)) {
      # Read the uploaded file
      df <- read.xlsx(input$dataFile$datapath)
      
      # Return the data frame to display as a table
      df
    }
  })
  
  # Generate a dendrogram from the uploaded data file
  output$moduleDendrogram <- renderPlot({
    # Check if a file is uploaded
    if (!is.null(input$dataFile)) {
      # Read the uploaded file
      df <- read.xlsx(input$dataFile$datapath)
      
      # Generate dendrogram
      dendro <- hclust(dist(df))
      
      # Plot the dendrogram
      plot(dendro, main = "Gene Co-expression Dendrogram")
    }
  })
  
  # Generate a module network plot from the uploaded data file
  output$moduleNetwork <- renderPlot({
    # Check if a file is uploaded
    if (!is.null(input$dataFile)) {
      # Read the uploaded file
      df <- read.xlsx(input$dataFile$datapath)
      
      # Perform network analysis and generate the module network plot
      
      # Convert the data frame to a numeric matrix
      datExpr <- as.matrix(df)
      
      # Create a signed adjacency matrix using the correlation-based approach
      adjacencyMatrix <- adjacency(datExpr, type = "signed", power = 6)
      
      # Perform network analysis using the adjacency matrix
      network <- blockwiseModules(adjacencyMatrix, power = 6)
      
      # Extract module assignments
      moduleColors <- network$colors
      
      # Create an igraph object from the adjacency matrix
      g <- graph_from_adjacency_matrix(adjacencyMatrix, mode = "undirected", weighted = TRUE)
      
      # Set the vertex color based on the module assignments
      V(g)$color <- moduleColors
      
      # Plot the module network
      plot(g, layout = layout_with_fr, vertex.size = 10,
           vertex.label = NA, vertex.label.cex = 0.8,
           main = "Module Network")
    }
  })
  
  # # Generate module assignments for the genes
  # output$moduleTable <- renderDataTable({
  #   # Check if a file is uploaded
  #   if (!is.null(input$dataFile)) {
  #     # Read the uploaded file
  #     df <- read.xlsx(input$dataFile$datapath)
  # 
  #     # Perform module detection using WGCNA
  #     exprsData <- as.data.frame(t(df))
  #     datExpr <- as.data.frame(exprsData)  # Convert to data frame
  #     geneTree <- hclust(dist(t(df)))
  #     modules <- cutreeDynamic(geneTree, method = "hybrid", deepSplit = 2)
  # 
  #     # Add module assignments to the data frame
  #     datExpr$Module <- modules
  # 
  #     # Return the data frame with module assignments
  #     datExpr
  #   }
  # })
}

# Run the application
shinyApp(ui, server)

