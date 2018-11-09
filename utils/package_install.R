cran_paks<-c("shiny","shinydashboard","plotly",
        "DT","ggplot2","plyr","dplyr",
        "shinytoastr","shinyWidgets","shinydashboardPlus",
        "shinycssloaders","shinytoastr","shinyjs","DBI",
        "reshape2","viridis", "BiocManager",
        "tm","data.table","RSQLite",
        "shinyBS","purrr", "devtools", "tidyverse")


github_paks<-c("nik01010/dashboardthemes")

bioconductor_paks<-c("GenomicFeatures")

already_installed<-rownames(installed.packages())

for (pak in cran_paks){
  if(pak in alread_installed){
	next
  } else {
  install.packages(pak, dependencies = T, repos = 'https://cloud.r-project.org/')
  }
}

for(pak in github_paks){
  devtools::install_github(pak)
}

for(pak in bioconductor_paks){
  BiocManager::install(pak, version = "3.8")
}
