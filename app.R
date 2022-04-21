library(shiny)
library(readxl)
library(DT)
library(bslib)
# library(dplyr)

# This sets the application-scoped cache to be a disk
# cache that can be shared among multiple concurrent R processes, and is
# deleted when the system reboots.
shinyOptions(cache = cachem::cache_disk(file.path(dirname(tempdir(
	
)), "WGCNA-cache")))

path_to_file <-  here::here("RefEx_expression_EST10_human.tsv")

# Load the data from the file
Database <- read.delim(path_to_file, row.names = 1)

# Define UI for application
ui <- fluidPage(# Theme
	theme = bs_theme(bootswatch = "simplex"),
	
	tags$html(
		tags$head(
			tags$meta(charset = "utf-8"),
			tags$meta(name = "description", content = "A Shiny App for WGCNA Analysis"),
			tags$meta(name = "robots", content = "noindex"),
			# Prevents app from showing up in search results
			# Favicon
			tags$link(
				rel = "shortcut icon",
				href = here::here("www", "favicon.ico"),
				type = "image/x-icon"
			)
		),
		
		tags$body(
			# Application title
			titlePanel("WGCNA Analysis"),
			h1("About"),
			div(
				p(
					"This is an online application that lets users perform",
					a("weighted gene co-expression network analysis", href = "https://en.wikipedia.org/wiki/Weighted_correlation_network_analysis"),
					"on their data."
				),
				p("To know more, take a look at the example below."),
				p(
					"It uses the R package ",
					a("WGCNA.", href = "https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/"),
					"A choice has been provided to the user in some cases. All other parameters have been left to their recommended defaults. "
				),
				p(
					"This is a ",
					a("Shiny Application.", href = "https://shiny.rstudio.com/")
				)
			),
			# The data table
			h1("Dataset"),
			DTOutput("GeneTable"),
			
			WGCNAShinyUI("1")
		),
		lang = "en"
	))

# Define server logic
server <- function(input, output) {
	output$GeneTable <- renderDT(
		Database,
		extensions = c("FixedColumns", "Buttons"),
		rownames = FALSE,
		options = list(
			searching = TRUE,
			searchHighlight = TRUE,
			autoWidth = TRUE,
			paging = TRUE,
			lengthMenu = matrix(
				c(10, 25, 50,-1, 10, 25, 50, "All"),
				ncol = 4,
				nrow = 2,
				byrow = TRUE
			),
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
			# buttons = c('copy', 'csv', 'excel', 'pdf', 'print'),
			buttons = c('copy', 'csv', 'excel'),
			dom = 'lfrtiBp'
		),
		server = TRUE
	)
	WGCNAShiny("1")
}

# Run the application
runApp(shinyApp(ui = ui, server = server),
	   host = "0.0.0.0",
	   port = 4330)
