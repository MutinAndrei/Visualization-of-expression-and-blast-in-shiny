This app helps you create a local blast with the right data and also allows you to track, visually, the level of expression.

Requirements for BLAST(https://github.com/ScientistJake/Shiny_BLAST):
Unix like environment
Current BLAST executables (https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
The following R packages:
library(DT)
library(thematic)
library(bslib)
library(XML)
library(plyr)
library(dplyr)
library(rBLAST)
library(ggplot2)
library(tidyr)
library(xlsx)
library(shinyalert)
library(Biostrings)
library(shinycssloaders)
library(ggtext)
library(plotly)
library(glue)
library(rclipboard)
library(BiocManager)
shinythemes (optional)

Note!
In order to run a local blast (https://github.com/ScientistJake/Shiny_BLAST/blob/master/README.md) you need to create your database. The second page allows you to create a table with transcript sequences, also made from your database