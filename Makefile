all: Stuff/pcounts.rds Pre_Implement.html

#Simulated population
Stuff/pcounts.rds: R/Simulated_Population.R R/Framework_functions.R
	Rscript R/Simulated_Population.R


#Preimplementation 
Pre_Implement.html: Pre_Implement.Rmd R/Framework_functions.R Stuff/pcounts.rds Stuff/removal_data.rds Stuff/removal_data_sub.rds Stuff/San_Nicolas_projected.shp
	Rscript -e "rmarkdown::render('Pre_Implement.Rmd')"

#Implementation


#Decision Support


#Shiny App

