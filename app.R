library(shiny)
library(bslib)
library(DT)
# library(dplyr)

# This sets the application-scoped cache to be a disk
# cache that can be shared among multiple concurrent R processes, and is
# deleted when the system reboots.
shinyOptions(cache = cachem::cache_disk(file.path(dirname(tempdir()), "WGCNA-cache")))

light <- bs_theme(bootswatch = "simplex")
dark <- bs_theme(
			bootswatch = "simplex",
			bg = "#292929",
			fg = "#999999",
			primary = "crimson",
			secondary = "coral"
		)
# dark <- bs_theme(bootswatch = "simplex", bg = "#204020", fg = "royalblue", primary = "turquoise", secondary = "teal")

# Define UI for application
ui <- fluidPage(
	# Theme
	theme = light,
	title = "WGCNA Analysis Online",
	
	tags$html(
		tags$head(
			tags$meta(charset = "utf-8"),
			tags$meta(name = "description", content = "A Free and Open Source app for online WGCNA Analysis"),
			# Favicon
			tags$link(
				rel = "shortcut icon",
				href = "favicon.png",
				type = "image/x-icon"
			)
		),
		
		tags$body(
			fluidRow(
				column(width = 10,
					h1("WGCNA Analysis Online")
				),
				column(width = 2,
					shinyWidgets::switchInput(
						inputId = "dark_mode",
						value = FALSE,
						onLabel = "Dark",
						offLabel = "Light",
						inline = TRUE
					)
				),
			),

			h2("About"),
			div(
				p(
					"This is an online application that lets users perform",
					a("weighted gene co-expression network analysis", href = "https://en.wikipedia.org/wiki/Weighted_correlation_network_analysis"),
					"on their data."
				),
				p(
					"To know more, take a look at the example below. Sample data was obtained from ",
					a("RefEx.", href = "http://refex.dbcls.jp/download.php?lang=en")
				),
				p(
					"It uses the R package ",
					a("WGCNA.", href = "https://horvath.genetics.ucla.edu/html/CoexpressionNetwork/Rpackages/WGCNA/"),
					"A choice has been provided to the user in some cases. All other parameters have been left to their recommended defaults. "
				),
				p(
					"This project is open-source. The source code is available on ",
					a("GitHub.", href = "https://github.com/Sciwhylab/WGCNA-Online")
				),
				pre(
					"ℹ This site is best viewed on large-screen devices.",
					class="alert alert-info"
				)
			),
			
			DataManagerUI("1"),
			# The data table
			h2("Data Viewer"),
			DTOutput("GeneTable"),
			WGCNAShinyUI("1")
		),
		lang = "en"
	)
)

# Define server logic
server <- function(input, output, session) {
	observe(session$setCurrentTheme(
		if (input$dark_mode)
			dark
		else
			light)
		)
	Database <- DataManager("1")
	WGCNAShiny("1", Database)
	output$GeneTable <- renderDT(
		Database(),
		extensions = c("FixedColumns", "Buttons"),
		rownames = FALSE,
		options = list(
			searching = TRUE,
			searchHighlight = TRUE,
			autoWidth = TRUE,
			paging = TRUE,
			lengthMenu = matrix(
				c(10, 25, 50, -1, 10, 25, 50, "All"),
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
}

# Run the application
## Local Only
# runApp(shinyApp(ui = ui, server = server),
# 	   host = "0.0.0.0",
# 	   port = 4330)
## Universal
shinyApp(ui = ui, server = server)
