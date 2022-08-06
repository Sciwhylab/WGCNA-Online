library(shiny)

# Setting up the env
options(stringsAsFactors = FALSE)

options(shiny.maxRequestSize = 50 * 1024 ^ 2) # Allow files uploads up to 50 MB
options(error=recover)

DataManagerUI <- function(id) {
	ns <- NS(id)
	
	tagList(
		h2("Data Source"),
		# textOutput(ns("Debug")), # Debugging only
		HTML(
			"<p>Sample data is <a href='https://figshare.com/articles/dataset/RefEx_expression_EST10_human_tsv_zip/4028625'>Processed expression data of 10 major tissues for EST human</a>.</p>"
		),
		actionButton(ns("useSample"), label = "Use Sample Data"),
		
		fileInput(
			ns("uploaded"),
			label = "Upload a file here",
			multiple = FALSE,
			accept = c(".csv", ".tsv", ".xlsx", ".xls")
		),
		checkboxGroupInput(
			ns("uploadOptions"),
			label = "Options for the Uploaded File",
			choices = c(
				"Contains Column Names" = "col_names",
				"Contains Row Names" = "row_names",
				"Interpret -1 as NA" = "m1toNA",
				"Normalization with log2" = "norm_log2"
			)
		),
		pre(
			"Different datasets require different treatment. It is recommended that you process your data before uploading it. After all, you know your data the best.",
			class="alert alert-warning"
		),
		actionButton(ns("useUploaded"), label = "Use My Data")
	)
}

DataManager <- function(id) {
	moduleServer(id,
				 function(input, output, session) {
				 	# Load the data from the sample file
				 	# File downloaded from http://refex.dbcls.jp/download.php?lang=en
				 	path_to_file <- here::here("sample_datasets/RefEx_expression_EST10_human.tsv")
				 	sampleData <- read.delim(path_to_file, row.names = 1)
				 	# Replace -1 by NA
				 	sampleData[sampleData == '-1'] <- NA
				 	# dim(sampleData)
				 	# names(sampleData)
				 	# View(sampleData)
				 	
				 	sampleData <- as.data.frame(lapply(sampleData, as.numeric))
				 	sampleData <- na.omit(sampleData)
				 	# Normalization with log2
				 	sampleData <- log2(sampleData)
				 	
				 	reactives <- reactiveValues(uploaded_data = NULL)
				 	# uploaded_data <- NULL
				 	# head(datExpr1[1:5, 1:5]) # samples in row, genes in column
				 	
				 	bindEvent(observe({
				 		reactives$uploaded_data <- NULL
				 	}),
				 	input$useSample)
				 	
				 	bindEvent(observe({
				 		file <- isolate(input$uploaded)
				 		ext <- tools::file_ext(file$datapath)
				 		
				 		req(file)
				 		uploadOptions <- isolate(input$uploadOptions)
				 		col_names <- "col_names" %in% uploadOptions
				 		if (ext == "csv")
				 		{
				 			reactives$uploaded_data <-
				 				readr::read_csv(file$datapath, col_names = col_names)
				 		} else if (ext == "tsv") {
				 			reactives$uploaded_data <-
				 				readr::read_tsv(file$datapath, col_names = col_names)
				 		} else if (ext == "xls") {
				 			reactives$uploaded_data <-
				 				readxl::read_xls(file$datapath, col_names = col_names)
				 		} else if (ext == "xlsx") {
				 			
				 		} else{
				 			validate(need(FALSE, "File format is not supported."))
				 		}
				 		
				 		if ("row_names" %in% uploadOptions)
				 			reactives$uploaded_data <- reactives$uploaded_data[,-1]
				 		
				 		if ("m1toNA" %in% uploadOptions)
				 			# Replace -1 by NA
				 			reactives$uploaded_data[reactives$uploaded_data == '-1'] <-
				 			NA
				 		
				 		reactives$uploaded_data <-
				 			lapply(reactives$uploaded_data, as.numeric) %>% 
				 			as.data.frame()
				 		
				 		# message(all.equal(sampleData,isolate(reactives$uploaded_data)))
				 		tryCatch({
				 			gsg = goodSamplesGenes(reactives$uploaded_data, verbose = 3);
				 			gsg$allOK
				 			if (!gsg$allOK)
				 			{
				 				# Optionally, print the gene and sample names that were removed:
				 				if (sum(!gsg$goodGenes)>0)
				 					printFlush(paste("Removing genes:", paste(names(reactives$uploaded_data)[!gsg$goodGenes], collapse = ", ")));
				 				if (sum(!gsg$goodSamples)>0)
				 					printFlush(paste("Removing samples:", paste(rownames(reactives$uploaded_data)[!gsg$goodSamples], collapse = ", ")));
				 				# Remove the offending genes and samples from the data:
				 				reactives$uploaded_data = reactives$uploaded_data[gsg$goodSamples, gsg$goodGenes]
				 			}
				 		},
				 		error = function(e) {
				 			message("Warning: Could not get good genes. Proceeding with NA removal.")
				 			reactives$uploaded_data <- reactives$uploaded_data %>% na.omit()
				 		})
				 		# message("Reached this")
				 		
				 		if ("norm_log2" %in% uploadOptions)
				 			# Normalization with log2
				 			reactives$uploaded_data <- log2(reactives$uploaded_data)
				 		
				 		output$Debug <- renderText({
				 			file$datapath
				 		})
				 		
				 	}),
				 	input$useUploaded)
				 	
				 	datExpr1 <- reactive({
				 		if (is.null(reactives$uploaded_data))
				 			sampleData
				 		else
				 			reactives$uploaded_data
				 	})
				 	
				 	return(datExpr1)
				 })
}
