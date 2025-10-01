
library(shiny)
library(tidyverse)

denver_tree_sf <- read_rds("../data_raw/denver_tree/denver_tree_clean.rds")


# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Denver Trees Exploritory Analyis Application"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
          dateRangeInput("data_range", "Select a Tree Inventrory Range?", 
                         start = min(denver_tree_sf$inventory_date),
                         end = max(denver_tree_sf$inventory_date),
                         min = min(denver_tree_sf$inventory_date),
                         max = max(denver_tree_sf$inventory_date)
                         ),
          varSelectInput("sci_species", )
        ),

        # Show a plot of the generated distribution
        mainPanel(
           plotOutput("map_denver")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$distPlot <- renderPlot({

    })
}

# Run the application 
shinyApp(ui = ui, server = server)
