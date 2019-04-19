library(shiny)
library(knitr)
library(raster)
library(sf)
library(ggspatial)
library(tidyverse)
library(R.utils)
library(kableExtra)
source("R/Framework_functions.R")

##DEFINE THE USER INTERFACE################################################################################################
ui<-fluidPage(
	#code to make a horizontal line in a sidebar
	tags$head(
		tags$style(HTML("hr {border-top: 1px solid #000000;}")),
		tags$style(HTML(".secondary .widget {margin-bottom: 0 !important;}"))
	),
	titlePanel("Monitoring of implementation"),
	sidebarLayout(
		sidebarPanel(#INPUT THE REGION, detector and detections FILEs ---------------------------------------------------
				   h3("Input data"),
				   fileInput(inputId="boundary", label="Region file", accept=c(".gpkg")),
				   fileInput(inputId="traps", label="Trap locations",   accept=c("text/csv", ".csv")),
				   fileInput(inputId="removals", label="Trapping histories", accept=c("text/csv", ".csv")),
				   fileInput(inputId="detectors", label="Passive detector locations",   accept=c("text/csv", ".csv")),
				   fileInput(inputId="detections", label="Detection histories", accept=c("text/csv", ".csv")),
				   hr(),
				   #GRAPH REMOVALS AND SIGN SURVEYS
				   h3("Plot data"),
				   actionButton("RemGraph", "Removal Graph"),
				   #SELECT APPROPRIATE MODEL --------------------------------------------------------------------------
				   h3("Inference"),
				   radioButtons(inputId="Data", label="Data Types", inline=FALSE,
				   		   choices=c("Removal Only"="remov", "Removal+Detections"="removsign", "Detections Only"="sign")),
				   #MCMC options and trigger fitting of model ---------------------------------------------------------
				   splitLayout(
				   	numericInput(inputId="burnin", label="burnin", min=500, max=10000, value=500, step=100),
				   	numericInput(inputId="samples", label="samples", min=500, max=10000, value=1000, step=100),
				   	numericInput(inputId="chains", label="chains", min=1, max=10, value=2, step=1)
				   ),
				   actionButton("Fit_mod", "Fit model"),
				   h3("Results"),
				   actionButton("PlotRemoval", "Plot removal"),
				   actionButton("PlotPeradicated", "Calculate Pr(eradication)")
		),
		mainPanel(
			plotOutput(outputId = "map"),
			verbatimTextOutput(outputId = "designstats")
		)
	)	
)


##DEFINE THE SERVER#####################################################################################################
server<-function(input, output){ 
	#reactive import of region file
	region <- reactive({ read_sf(input$boundary$datapath) })
	
	#reactive import of trap locations
	traps <- reactive({ read_sf(input$traps$datapath) })
	
	#reactive import of trap removals
	removals <- reactive({ read_csv(input$removals$datapath) })
	
	#reactive import of detector locations
	detectors <- reactive({ read_csv(input$detectors$datapath) })
	
	#reactive import of detection histories
	detections <- reactive({ read_csv(input$detections$datapath) })
	#----Output graph of removal progress
        #code for graph of cumulative removals, and detections of sign over time.
	
	#----Output map (animated perhaps) of trapping and detections
	   #animated map of trapping and detection
	
	#----model fit and results
	   #graph of estimated abundance or occupancy over time.
	   #
}

shinyApp(ui=ui, server=server)