library(shiny)
library(shinyWidgets)
library(tidyverse)
library(leaflet)
library(sf)
library(tigris)
library(DT)

denver_tree_sf <- read_rds("../data_raw/denver_tree/denver_tree_clean.rds")

denver_county <- counties(state = "CO", cb = TRUE) %>%
  filter(NAME == "Denver") %>%
  st_transform(crs = 4326) # adjust to local projection for final sample

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Denver Trees Exploritory Analyis Application"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
          h5("Tree Inputs"),
          dateRangeInput("data_range", "Select a Date Range?", 
                         start = min(denver_tree_sf$inventory_date),
                         end = max(denver_tree_sf$inventory_date),
                         min = min(denver_tree_sf$inventory_date),
                         max = max(denver_tree_sf$inventory_date)
                         ),
          sliderTextInput("diameter", 
                          "Diameter of the Tree",
                          choices = c("0 to 6", "6 to 12", "12 to 18", 
                                      "18 to 24", "24 to 30", "36 to 42", 
                                      "42 to 48", "48 +"),
                          selected = c("0 to 6", "48 +"),
                          grid = TRUE),
          hr(),
          h5("Map Inputs"),
          switchInput("toggle_map_view", 
                      label = "Map View",
                      value = FALSE,
                      onLabel = "Satellite",
                      offLabel = "Topographic",
                      size = "mini"),
          switchInput("toggle_map_obs", 
                      label = "Map Obs.",
                      value = FALSE,
                      onLabel = "Individual",
                      offLabel = "Cluster",
                      size = "mini")
          ),
        # Show a plot of the generated distribution
        mainPanel(
          leafletOutput("map_denver")
           )
        ),
    downloadButton("download_obs", "Download All Selected Observations"),
    fluidRow(
      column(12,
             DTOutput("tree_table")
             )
      ),
    conditionalPanel(
      condition = "input.tree_table_rows_selected.length > 0",
      hr(),
      h3("Selected Trees - Click to Expand Details"),
      uiOutput("selected_accordions")
    )
    )

# Define server logic required to draw a histogram
server <- function(input, output, session) {  # Add 'session' parameter
  
  filtered_trees <- reactive({
    denver_tree_sf %>%
      filter(inventory_date >= input$data_range[1],
             inventory_date <= input$data_range[2]) %>%
      filter(diameter >= input$diameter[1],
             diameter <= input$diameter[2])
  })
  
  grouped_data <- reactive({
    filtered_trees() %>%
      group_by(species_sci, 
               species_common, 
               dc_tree, 
               leaf_persistence, 
               growth_habit) %>% 
      summarise(Obs. = n(), .groups = "drop")
    })
  
  output$tree_table <- renderDT({
    datatable(
      grouped_data(),
      selection = list(mode = "multiple", selected = NULL),
      options = list(pageLength = 10),
      rownames = FALSE
      )
    })
  
  selected_trees <- reactive({
    if (length(input$tree_table_rows_selected) > 0) {
      selected_rows <- grouped_data()[input$tree_table_rows_selected, ]
      filtered_trees() %>%
        filter(species_sci %in% selected_rows$species_sci)
    } else {
      filtered_trees() %>% 
        filter(FALSE)  # empty sf
    }
  })
  
  output$map_denver <- renderLeaflet({
    leaflet() %>%
      addProviderTiles(providers$Esri.WorldTopoMap) %>%
      addPolygons(
        data = denver_county,
        fillOpacity = 0,
        color = "red",
        weight = 2,
        opacity = 0.8
      )
  })
  # Toggle Map View
  observeEvent(input$toggle_map_view, {
    leafletProxy("map_denver")%>% 
      clearTiles() -> proxy
    
    if(input$toggle_map_view) {  # TRUE = Topographic
      proxy %>% 
        addProviderTiles(providers$Esri.WorldImagery)
      } else {  # FALSE = Satellite
        proxy %>% 
          addProviderTiles(providers$Esri.WorldTopoMap) 
      }
    })
  
  # Observe: update map when selection changes
  observeEvent(list(input$toggle_map_obs, input$tree_table_rows_selected), {
    proxy <- leafletProxy("map_denver")
    
    # Clear old layers
    proxy %>%
      clearMarkers() %>%
      clearMarkerClusters()
    
    trees <- selected_trees()
    
    if (nrow(trees) == 0) return()  # nothing selected → keep map empty
    
    if (input$toggle_map_obs) {
      # Individual points
      proxy %>%
        addCircleMarkers(
          data = trees,
          radius = 2,
          color = "red",
          stroke = FALSE,
          fillOpacity = 1
        )
    } else {
      # Clustered points
      proxy %>%
        addCircleMarkers(
          data = trees,
          radius = 3,
          color = "red",
          stroke = FALSE,
          fillOpacity = 1,
          clusterOptions = markerClusterOptions()
        )
    }
  })
  
  
  output$download_obs <- downloadHandler(
    filename = function() {
      paste0("denver_trees_selected_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(input$tree_table_rows_selected)
      
      # Get selected species
      selected_rows <- grouped_data()[input$tree_table_rows_selected, ]
      
      # Filter and prepare data
      download_data <- filtered_trees() %>%
        filter(species_sci %in% selected_rows$species_sci) %>%
        st_drop_geometry() %>%  # Remove geometry
        select(-c(globalid, objectid))  # Optional: remove ID columns
      
      # Write to CSV
      write.csv(download_data, file, row.names = FALSE)
    }
  )
  
  # output$download_obs <- downloadHandler(
  #   filename = function() {
  #     paste0("denver_trees_", Sys.Date(), ".geojson")
  #   },
  #   content = function(file) {
  #     req(input$tree_table_rows_selected)
  #     
  #     selected_rows <- grouped_data()[input$tree_table_rows_selected, ]
  #     
  #     download_data <- filtered_trees() %>%
  #       filter(species_sci %in% selected_rows$species_sci)
  #     
  #     st_write(download_data, file, driver = "GeoJSON", delete_dsn = TRUE)
  #   }
  # )
  
}

# Run the application 
shinyApp(ui = ui, server = server)
