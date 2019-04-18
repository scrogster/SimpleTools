all: Stuff/pcounts.rds Pre_Implement.html Implementation.html

#Simulated population
Stuff/pcounts.rds: R/Simulated_Population.R R/Framework_functions.R
	Rscript R/Simulated_Population.R


#Preimplementation 
Pre_Implement.html: Pre_Implement.Rmd R/Framework_functions.R Stuff/pcounts.rds Stuff/removal_data.rds Stuff/removal_data_sub.rds Stuff/San_Nicolas_projected.shp
	Rscript -e "rmarkdown::render('Pre_Implement.Rmd')"

#Implementation
Implementation.html: Implementation.Rmd R/Framework_functions.R Pre_Implement.html
	Rscript -e "rmarkdown::render('Implementation.Rmd')"

#Decision Support


#Shiny App

