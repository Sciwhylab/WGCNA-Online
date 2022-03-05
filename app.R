
library(shiny)
library(readxl)
library(DT)
library(bslib)
# library(dplyr)

path_to_file <-  here::here("RefEx_expression_EST10_human.tsv")

# Load the data from the file
Database <- read.delim(path_to_file, row.names = 1)

# Define UI for application
ui <- fluidPage(
	# Theme
	theme = bs_theme(
		bootswatch = "minty"
	), 
	# # Favicon
	# tags$head(tags$link(rel="shortcut icon", href=here::here("www", "favicon.ico"), type="image/x-icon")),
	# Application title
	titlePanel("RefEx_expression_EST10_human"), 
	# The data table
	DTOutput("GeneTable")
)

# Define server logic
server <- function(input, output) {
	
	output$GeneTable <- renderDT(Database, 
								 extensions = c("FixedColumns", "Buttons"), 
								 rownames = FALSE, 
								 options = list(
								 	searching = TRUE,
								 	searchHighlight = TRUE,
								 	autoWidth = TRUE,
								 	paging = TRUE,
								 	lengthMenu = matrix(c(10, 25, 50, -1, 10, 25, 50, "All"), ncol = 4, nrow = 2, byrow = TRUE),
								 	pagingType = "full_numbers",
								 	scroller = TRUE,
								 	scrollX = TRUE,
								 	scrollY = "900px",
								 	scrollCollapse = TRUE,
								 	fixedHeader = TRUE,
								 	class = 'cell-border stripe',
								 	# fixedColumns = list( # unnecessary because the 
								 	# 	leftColumns = 1,
								 	# 	heightMatch = 'none'
								 	# ),
								 	buttons = c('copy', 'csv', 'excel', 'pdf', 'print'), 
								 	dom = 'lfrtiBp'
								 ), 
								 server = TRUE
	)
}

# Run the application 
runApp(shinyApp(ui = ui, server = server), host = "0.0.0.0", port = 4330)
