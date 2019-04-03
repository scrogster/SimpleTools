SimpleTools - Decision support tools to support pest eradication
projects
================

**File descriptions**:

*Simulated Population.R* – R code used to generate test data from a
simulated pest eradication. The simulated data is saved in Stuff/. for
use in *Pre\_Implement.Rmd* or *Implementation.Rmd*.

*Framework\_functions.r* - Functions and models used in all \*Rmd files
or *Simulated Population.r* are found here.

*Pre\_Implement.Rmd* – R markdown document that undertakes an initial
assessment of the pest population to be eradicated, primarily to
estimate initial abundance. Input monitoring data and region shapefile
are required. Input monitoring data can be generated with *Simulated
Population.r*

*Implementation\_Rmd* – R markdown document that takes data on the
sequential removal of individuals from the pest population and uses it
to estimate the residual population remaining. Can also make use of
additional, non-removal monitoring data. Initially this is confined to
camera monitoring or sand plot (sign) monitoring. Required data can be
generated using *Simulated Population.r*.

*Decision\_Support.Rmd* - R markdown document that quantifies the
probability of eradication success from monitoring that detects no
individuls of the pest species. Designed to be used at the point when it
is thought that eradication may have been achieved. Uses a Bayesian
framework and detection parameters to estimate the amount of monitoring
required with zero detections to achieve eradication through a “stopping
rule”. Currently the stopping rule is a threshold on the type I error
rate (probability of falsely declaring eradication successful)

## Prerequisites

The scripts require packages *tidyverse*, *raster*, *spatstat*, *sf* and
*maptools*.
