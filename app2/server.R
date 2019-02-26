

server<-function(input, output, session){
  
  source("../modules/gene_select.R")
  source("../modules/geneviz.R")
  source("../modules/tissue_select.R")
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
  
  output$login_menu<-renderUI({
    if(login_status$login){
      fluidRow(
        actionButton(inputId = "login_button", label = "Logout", icon=icon("sign-out"))
      )
    } else {
      fluidRow(
        textInput(inputId = "username", label = "Username"),
        passwordInput(inputId = "password", label = "Password"),
        actionButton(inputId = "login_button", label = "Login", icon=icon("sign-in"))
      )
    }
  })
  
  observeEvent(input$login_button, {
    if(login_status$login){
      login_status$login<-F
      toastr_success("Logout successful")
    } else {
      if(input$username=="User" & input$password=="pass"){
        toastr_success("Login successful")
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
  
  #this returns a data frame of gene
  gene_selection<-callModule(module=gene_select_server, id="gene_select", conn=gtex, login=login_status)
  
  #TODO return tissues and samples
  tissue_selection<-callModule(module=tissue_select, id="tissue_select", conn=gtex, login=login_status)
  
  callModule(module=median_heatmap, id="median_heatmap", conn=gtex, login=login_status, genesymbol=gene_selection)
  
  clicked_gene<-callModule(gene_expression, id="gene_expression", conn=gtex, login=login_status, 
                           genes=gene_selection, tissues_samples=tissue_selection)
  
  #TODO add the gene title in the module
  #callModule(module = genevis, id = "gtex_plot", 
  #           tissues_samples=tissue_selection,
  #           conn=gtex, gene_id=clicked_gene)
  
  session$onSessionEnded(
    function(){
      dbDisconnect(gtex)
    }
  )  
}