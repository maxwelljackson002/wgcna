library(shiny)
library(openxlsx)
library(WGCNA)
library(igraph)
library(BiocManager)
library(tidyverse)
library(shinyjs)
library(gplots)

options(repos = BiocManager::repositories())

# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("Weighted Gene Co-Expression Network Analysis"),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      fileInput("dataFile", "Upload Gene Expression Data", accept = c(".csv", ".xlsx")),
      sliderInput("numPointsDendrogram", "Number of Points (Dendrogram):", min = 1, max = 100, value = 30),
      sliderInput("numPointsNetwork", "Number of Points (Network):", min = 1, max = 100, value = 30),
      sliderInput("numPointsHeatmap", "Number of Points (Heatmap):", min = 1, max = 100, value = 30),
      actionButton("resetButton", "Reset")
    ),
    
    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(
        tabPanel("Module Assignment", dataTableOutput("moduleTable")),
        tabPanel("Module Dendrogram", plotOutput("moduleDendrogram")),
        tabPanel("Module Network", plotOutput("moduleNetwork")),
        tabPanel("Module Heatmap", plotOutput("heatmapPlot"))
      )
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
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
  
  # Generate a dendrogram from the uploaded data file with variable number of points
  output$moduleDendrogram <- renderPlot({
    # Check if a file is uploaded
    if (!is.null(input$dataFile)) {
      # Read the uploaded file
      df <- read.xlsx(input$dataFile$datapath)

      # Retrieve the number of points to display from the sliderInput
      num_points <- input$numPointsDendrogram

      # Generate dendrogram with the specified number of points
      dendrogram <- hclust(dist(df))

      # Cut the dendrogram into the desired number of clusters based on the number of points selected by the user
      dendrogram_cut <- cutree(dendrogram, k = num_points)

      # Plot the dendrogram with limited points
      plot(dendrogram, main = "Gene Co-expression Dendrogram with Limited Points")
      rect.hclust(dendrogram, k = num_points, border = 2)
    }
  })
  
  # Generate a module network plot from the uploaded data file
  output$moduleNetwork <- renderPlot({
    # Check if a file is uploaded
    if (!is.null(input$dataFile)) {
      # Read the uploaded file
      df <- read.xlsx(input$dataFile$datapath)
      
      # Retrieve the number of points to display from the sliderInput
      num_points <- input$numPointsNetwork
      
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
  
  # Generate a heatmap from the uploaded data file
  output$heatmapPlot <- renderPlot({
    # Check if a file is uploaded
    if (!is.null(input$dataFile)) {
      # Read the uploaded file
      df <- read.xlsx(input$dataFile$datapath)
      
      # Convert non-numeric values to NA in the entire data frame
      df[] <- lapply(df, as.numeric)
      
      # Retrieve the number of points to display from the sliderInput
      num_points <- input$numPointsHeatmap
      
      # Generate the heatmap using the heatmap.2 function from gplots
      heatmap.2(as.matrix(df), Rowv = NA, Colv = NA, scale = "row",
                trace = "none", col = rev(heat.colors(256)))
    }
  })
}

# Reactive value to store the uploaded data
uploaded_data <- reactiveVal(NULL)

# Function to reset the app's state
observeEvent(input$resetButton, {
  uploaded_data(NULL)
  output$moduleDendrogram <- NULL
  output$moduleNetwork <- NULL
  output$heatmapPlot <- NULL
})

# Process the uploaded file and store the data
observeEvent(input$dataFile, {
  # Check if a file is uploaded
  if (!is.null(input$dataFile)) {
    # Read the uploaded file
    df <- read.xlsx(input$dataFile$datapath)
    
    # Store the data in the reactive value
    uploaded_data(df)
  }
})

options(shiny.maxRequestSize = 50*1024^2)  # Set maximum file size to 50MB

# Run the application
shinyApp(ui, server)
