

server<-function(input, output, session){
  
  source("../modules/gene_select.R")
  source("../modules/geneviz.R")
  source("../modules/sample_subject_filter.R")
  source("../utils/getdata.R")
  
  toastr_info(message = "loading dependencies please wait...")
  
  suppressPackageStartupMessages(library(shiny))
  suppressPackageStartupMessages(library(plotly))
  suppressPackageStartupMessages(library(DT))
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(plyr))
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(shinytoastr))
  suppressPackageStartupMessages(library(shinyWidgets))
  suppressPackageStartupMessages(library(shinycssloaders))
  suppressPackageStartupMessages(library(shinytoastr))
  suppressPackageStartupMessages(library(shinyjs))
  suppressPackageStartupMessages(library(DBI))
  suppressPackageStartupMessages(library(reshape2))
  suppressPackageStartupMessages(library(viridis))
  suppressPackageStartupMessages(library(tm))
  suppressPackageStartupMessages(library(data.table))
  suppressPackageStartupMessages(library(shinyBS))
  
  options(warn=-1)
  
  credentials<-c("alper", "pass")
  gtex<-dbConnect(drv = RPostgreSQL::PostgreSQL(), host="localhost", 
                  dbname="gtex",
                  user=credentials[1], password=credentials[2])
  
  login_status<-reactiveValues(login=F)
  
  observeEvent(input$login_button, {
    if(login_status$login){
      login_status$login<-F
      toastr_success("Logout successful")
      updateActionButton(session, "login_button", label="Login")
    } else {
      if(input$username=="User" & input$password=="pass"){
        toastr_success("Login successful")
        updateActionButton(session, "login_button", label="Logout")
        login_status$login<-T
      } else {
        toastr_error("Login failed")
        login_status$login<-F
      }
    }
  })
  
  output$is_loggedin<-renderText({
    login_status$login
  })
  
  
  gene_tissue_selection<-callModule(module=gene_select, id="gene_select", conn=gtex, login=login_status)
  
  output$tissue_select_check<-renderText({
    unlist(gene_tissue_selection)
  })
  
  #callModule(median_heatmap, gene_tissue_selection$genes_df)
  
  #TODO add the gene title in the module
  callModule(module = genevis, id = "gtex_plot", 
             tissues=tissues, samples=filtered_samples,
             conn=gtex, gene_id=clicked_gene)
  
  session$onSessionEnded(
    function(){
      dbDisconnect(gtex)
    }
  )  
}