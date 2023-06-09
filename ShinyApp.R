#' ####################################################################### #
#' PROJECT: [Biodiversity Simplification & Ecological Network Topology] 
#' CONTENTS: 
#'  - Network Visualisations in Shiny App
#'  DEPENDENCIES:
#'  - [Executed] 3 - Visualisation.R
#' AUTHOR: [Erik Kusch]
#' ####################################################################### #
library(shiny)
library(shinythemes)
library(shinyWidgets)
library(ggplot2)
library(sf)
library(cowplot)

## CRAN -------------------------------------------------------------------
install.load.package <- function(x) {
  if (!require(x, character.only = TRUE))
    install.packages(x, repos='http://cran.us.r-project.org')
  require(x, character.only = TRUE)
}

package_vec <- c(
  "devtools", # needed for non-cran packages further down
  "rgeos", # for loading shapefiles
  "tidyverse", # for data handling
  "rgbif", # for occurrence retrieval
  "pbapply", # for apply with progress bar
  "data.table", # for data handling
  "rnaturalearth", # for landmask in projection kriging
  "rnaturalearthdata", # for landmask in projection kriging
  "rredlist", # for IUCN risk retrieval
  "ConR", # for computation of IUCN risks
  "CoordinateCleaner", # for additional occurrence cleaning
  "igraph", # for graph operations
  "FD", # for gower distance of trait data
  "reshape2", # for making network matrices into plottable data frames via melt()
  "bipartite", # for bipartite network analyses
  "leaflet", # for html map products to investigate networks separately
  "leafpop", # for graph popups in leaflet output
  "cowplot", # for arranging of plots
  "gridExtra", # for table grobs as legends in plots
  "dplyr" # for data cleaning
)
sapply(package_vec, install.load.package)


# PREAMBLE ==================================================================
# source("3 - Visualisation.R")
# source("0 - Preamble.R")
load(file = "Shiny.RData")

# SHINY APP =================================================================
## User Interface -----------------------------------------------------------
ui <- fluidPage(theme = shinytheme("slate"),
                navbarPage(
                  # theme = "cerulean",  # <--- To use a theme, uncomment this
                  "Biodiversity Simplification & Ecological Network Topology",
                  tabPanel("Network Measures",
                           mainPanel(
                             # h1("Filtering Criteria"),
                             # fluidRow(
                             #   column(3,
                             #          hr(),
                             #          pickerInput("Holding", "Holding", choices = Locations_sp$Holding,
                             #          options = list(`actions-box` = TRUE),
                             #          multiple = TRUE, selected = NULL)
                             #   ),
                             #   column(3,
                             #          hr(),
                             #          pickerInput("Species", "Species", choices = unique(Data_df$Species),
                             #                      options = list(`actions-box` = TRUE),
                             #                      multiple = TRUE, selected = NULL)
                             #          ),
                             # column(3,
                             #        hr(),
                             #        pickerInput("Location", "Location", choices = Locations_sp$ID,
                             #                    options = list(`actions-box` = TRUE),
                             #                    multiple = TRUE, selected = NULL),
                             #        ),
                             # ),
                             h1("Map Selection"),
                             leafletOutput("map", height = 700),
                             # absolutePanel(top = 10, right = 10,
                             #               pickerInput("countries", label = "Select a Country:",
                             #                           choices = list("All countries",
                             #                                          unique(leaflet_df$locality)
                             #                                          ),
                             #                           options = list(
                             #                             `live-search` = TRUE)
                             #                           )
                             #               ),
                             plotOutput("Plot_Net"),
                             # uiOutput("plot.ui")
                          ), # mainPanel
                           sidebarPanel(
                             tags$h1("Network Data"),
                             textOutput("nObs"),
                             tags$h2(textOutput("Study")),
                             tableOutput("NetTable"),
                             # fluidRow(verbatimTextOutput("map_marker_click")), 
                             tags$h2("Modularity"),
                             plotOutput("Plot_Mod"),
                             br(),
                             tags$h2("Nestedness"),
                             plotOutput("Plot_Nest")
                             ), # sidebarPanel
                           
                  ), # Data Explorer Panel
                  tabPanel("Species Abbreviations", 
                           fluidRow(
                             column(width = 5,
                                    tags$h1("Animals"),
                                    tableOutput("Tab_Animals")
                             ),
                             column(width = 5,
                                    tags$h1("Plants"),
                                    tableOutput("Tab_Plants")
                             ),
                           )
                           ),
                  tabPanel("The Team", "This panel is intentionally left blank")
                ) # navbarPage
) # fluidPage

## Server Function ----------------------------------------------------------
server <- function(input, output) {
  ### Behind the Scenes Data Manipulation ####
  #### Filtering Data ----
  ## this shiny function subsets the raw data according to user-specific data filters, currently not implemented
  Data.Subset <- reactive({
    SubsetData <- leaflet_df
    return(SubsetData)
  })
  
  ## retaining only plots queried by user
  Plots.Subset <- reactive({
    SubsetPlots <- Plots_ls
    SubsetPlots <- SubsetPlots[which(names(SubsetPlots) %in% Data.Subset()$net.id)]
    return(SubsetPlots)
  })
  
  ### Output for Shiny ####
  #### Main Panel ----
  output$map <- renderLeaflet({
    
    plot_df <- Data.Subset()
    
    netpaths <- file.path(Dir.E.LeafletImages, paste0(plot_df$net.id, ".png"))
    label_ls <- lapply(paste( "<b> Study ID: </b>" , plot_df$study.id, "<br>",
                              "Network ID:" , plot_df$net.id, "<br>"
                              , "Weighted Nestedness", plot_df$Nestedness, "<br>"
                              , "Modularity", plot_df$Modularity, "<br>",
                              "Robustness to animal extinction", plot_df$RobustnessAnimal, "<br>"
                              , "Robustness to plant extinction", plot_df$RobustnessPlant),
                       htmltools::HTML)
    
    Map <- leaflet() %>%
      # Base groups
      addTiles() %>%
      # Overlay groups
      addCircleMarkers(data = plot_df, ~longitude, ~latitude,
                       group = "net.id", 
                       label = label_ls
      ) 
    # %>%
    #   addPopupImages(netpaths, group = "net.id", width = 600, maxWidth = 600, height = 300)
    # addPopupGraphs(Plots.Subset(), group = "net.id", width = 800, maxWidth = 800, height = 400)
    return(Map)
    
  })
  
  #### Reactive Clicks ----
  clickedMarker <- reactiveVal()
  clickedMap <- reactiveVal()
  
  ## on click of markers in leaflet map
  observeEvent(input$map_marker_click,{
    # subsetting data
    event <- input$map_marker_click
    extract_df <- Data.Subset()[Data.Subset()$latitude == event$lat & Data.Subset()$longitude == event$lng, ]

    # adding newly coloured circle marker    
    proxy <- leafletProxy("map")
    proxy <- proxy %>%
      addCircleMarkers(group = as.character(extract_df$net.id),
                        lng=extract_df$longitude,
                        lat=extract_df$latitude,
                        color = "red")
    # removing previous red marker if one is present
    if(!is.null(clickedMarker())){
      prox <- proxy %>%
       clearGroup(as.character(clickedMarker()))
    }
    # saving objects
    clickedMarker(extract_df$net.id)
    proxy
    })
  
  observeEvent(input$map_click,{
    proxy <- leafletProxy("map")
    # removing previous red marker if one is present
    if(!is.null(clickedMap())){
      prox <- proxy %>%
        clearGroup(as.character(clickedMap()))
    }
    clickedMap(clickedMarker())
    clickedMarker("Empty")
  })
  
  #### Side Panel ----
  output$nObs <- renderText(paste("Showing", nrow(Data.Subset()), "individual networks.")) # for reporting number of data points
  # this is the visualisation of the aggregate PI in the sample with corresponding colouring
  output$Plot_Mod <- renderPlot({
    plot_df <- Data.Subset()
    ggplot(plot_df, aes(x = Modularity)) +
      geom_histogram(bins = 50) + 
      geom_density(size = 2) +
      theme_bw() + labs(y = "Count")  + 
      labs(title = "Pre-Extinction Modularity")
  })
  
  output$Plot_Nest <- renderPlot({
    plot_df <- Data.Subset()
    ggplot(plot_df, aes(x = Nestedness)) +
      geom_histogram(bins = 50) + 
      geom_density(size = 2) +
      theme_bw() + labs(y = "Count")  + 
      labs(title = "Pre-Extinction Nestedness")
  })
  
  output$Plot_Net <- renderPlot({ggplot() + theme_void()})
  
  output$Tab_Animals <- renderTable(AnimalsLegend)
  output$Tab_Plants <- renderTable(PlantsLegend)
  output$NetTable <- renderTable(
    data.frame(Modularity = mean(Data.Subset()$Modularity),
               Nestedness = mean(Data.Subset()$Nestedness))
  )
  
  
  
  observeEvent(clickedMap(), {
    event <- clickedMap()
    print(event)

    if(event != "Empty"){
      extract_df <- Data.Subset()[Data.Subset()$net.id %in% event, ]

      output$Study <- renderText(paste(unique(extract_df$study.id), collapse = " & "))
      output$NetTable <- renderTable(
        data.frame(Modularity = extract_df$Modularity,
                   Nestedness = extract_df$Nestedness)
      )

      output$Plot_Nest <- renderPlot({
        plot_df <- Data.Subset()
        ggplot(plot_df, aes(x = Nestedness)) +
          geom_histogram(bins = 50) +
          geom_density(size = 2) +
          geom_vline(xintercept = extract_df$Nestedness, size = 1.5, col = "red") +
          theme_bw() + labs(y = "Count")  +
          labs(title = "Pre-Extinction Nestedness")
      })

      output$Plot_Mod <- renderPlot({
        plot_df <- Data.Subset()
        ggplot(plot_df, aes(x = Modularity)) +
          geom_histogram(bins = 50) +
          geom_density(size = 2) +
          geom_vline(xintercept = extract_df$Modularity, size = 1.5, col = "red") +
          theme_bw() + labs(y = "Count")  +
          labs(title = "Pre-Extinction Modularity")
      })

      output$Plot_Net <- renderPlot({
        plot_grid(plotlist = Plots.Subset()[names(Plots.Subset()) %in% extract_df$net.id], ncol = 1)
      }, height = nrow(extract_df)*400)

    }else{
      
      output$Study <- renderText("")
      output$NetTable <- renderTable(
        data.frame(Modularity = mean(Data.Subset()$Modularity),
                   Nestedness = mean(Data.Subset()$Nestedness))
      )

      output$Plot_Mod <- renderPlot({
        plot_df <- Data.Subset()
        ggplot(plot_df, aes(x = Modularity)) +
          geom_histogram(bins = 50) +
          geom_density(size = 2) +
          theme_bw() + labs(y = "Count")  +
          labs(title = "Pre-Extinction Modularity")
      })

      output$Plot_Nest <- renderPlot({
        plot_df <- Data.Subset()
        ggplot(plot_df, aes(x = Nestedness)) +
          geom_histogram(bins = 50) +
          geom_density(size = 2) +
          theme_bw() + labs(y = "Count")  +
          labs(title = "Pre-Extinction Nestedness")
      })

      output$Plot_Net <- renderPlot({ggplot() + theme_void()})
      
    }



    # output$plot.ui <- renderUI({plotOutput("Plot_Net", height = 300*nrow(extract_df))})

  })
  
  
  
} # server

## Shiny Launch -------------------------------------------------------------
shinyApp(ui = ui, server = server)