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
	titlePanel("Pre-implementation Monitoring"),
	sidebarLayout(
		sidebarPanel("Data input",
				   #INPUT THE REGION, detector and detections FILEs ---------------------------------------------------
				   fileInput(inputId="boundary", label="Region file", accept=c(".gpkg")),
				   fileInput(inputId="detectors", label="Detector locations",   accept=c("text/csv", ".csv")),
				   fileInput(inputId="detections", label="Detection histories", accept=c("text/csv", ".csv")),
				   numericInput(inputId="res", label="Grid cell size", min=0, max=NA, value=1000, step=100),
				   radioButtons(inputId="gridtype", label="Grid type", c("Square"=1, "Hex"=0)),
				   #SELECT APPROPRIATE MODEL --------------------------------------------------------------------------
				   radioButtons(inputId="Model", label="Choose Model", choices=c("N-mixture"="Nmix", "Royle-Nichols"="RN")),
				   #MCMC options and trigger fitting of model ---------------------------------------------------------
				   numericInput(inputId="burnin", label="burnin", min=500, max=10000, value=500, step=100),
				   numericInput(inputId="samples", label="samples", min=500, max=10000, value=1000, step=100),
				   actionButton("Fit_mod", "Fit model")
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
	
	#reactive import of detector locations
	detectors <- reactive({ read_csv(input$detectors$datapath) })
	
	#reactive import of detection histories
	detections <- reactive({ read_csv(input$detections$datapath) })
	
 #----Output and render map of study area.
	output$map<-renderPlot({
		#requirements and imports
		req(input$boundary)
		req(input$detectors)
		req(input$gridtype)
		reg<-region()
		#make detection device locations as sf
		detectors<-st_as_sf(detectors(), coords = c("x", "y"), crs=st_crs(reg))
		#nd is the number of cells with a detection device
		nd<-length(detectors)
		gtype=input$gridtype==1
		#full grid
		tessel<-st_make_grid(reg, cellsize=input$res, square=gtype)
		tessel_cents<-st_make_grid(reg, cellsize=input$res, square=gtype, what = "centers") 
		a<-st_union(reg) %>% st_intersects(tessel) %>% as_tibble() 
		tessel<-tessel[a$col.id]
		tessel_cents<-tessel_cents[a$col.id]
		
		p1<- ggplot() +
			geom_sf(data= reg, fill="green", alpha=0.3) +
			geom_sf(data=tessel, fill=NA) +
			geom_sf(data=detectors, col="black") +
			annotation_scale(location = "bl", width_hint = 0.4) +
			annotation_north_arrow(location = "bl", which_north = "true", 
							   pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
							   style = north_arrow_fancy_orienteering) +
			labs(x="Longitude",y="Latitude",title="Sample region")+
			coord_sf()+
			theme_bw()
		p1}, 
		execOnResize = TRUE)	
    
}

shinyApp(ui=ui, server=server)