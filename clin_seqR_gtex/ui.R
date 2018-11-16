###### Clin SeqR Shiny app ######

# this is a two file shinly application the other
# one is called server.R

# Alper Celik



suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(DT))
suppressPackageStartupMessages(library(shinytoastr))
suppressPackageStartupMessages(library(shinyWidgets))
suppressPackageStartupMessages(library(shinycssloaders))
suppressPackageStartupMessages(library(shinytoastr))
suppressPackageStartupMessages(library(shinyjs))
suppressPackageStartupMessages(library(DBI))
suppressPackageStartupMessages(library(RSQLite))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(shinyBS))
suppressPackageStartupMessages(library(shinythemes))


source("../modules/geneviz.R")
source("../modules/sample_subject_filter.R")

median_connect<-dbConnect(drv=RSQLite::SQLite(), dbname="../data/namedb.db")

ui<-navbarPage("ClinSeqR", 
               tabPanel("Explore GTEx Data", 
                        fluidRow(
                          sidebarLayout(
                            sidebarPanel(width=3,
                                         tags$h3("Select Gene(s)"),
                                         radioGroupButtons(inputId = "gtex_selection_type", 
                                                           label = "", 
                                                           choices = c("Select From List", "Enter Text", "Select Gene Panel"), 
                                                           justified = TRUE,
                                                           selected = "Enter Text",individual = F, direction = "vertical"),
                                         withSpinner(
                                           uiOutput("gene_select"), color = '#8D8D8D'
                                         ), 
                                         br(),
                                         pickerInput("selected_tissues", label = "Select Tissues",
                                                     choices = dbGetQuery(
                                                       median_connect, "PRAGMA table_info(median_expression)")[2][-c(1:2),], 
                                                     options=list(`max-options`=10, 
                                                                  `live-search`=T, size=10), 
                                                     multiple = TRUE), 
                                         bsModal("sample_filter_modal", title = "Filter GTEx Data", 
                                                 trigger = "filter_gtex_data", size = "large", 
                                                 filter_modal_ui("gtex_filter")),
                                         actionButton("filter_gtex_data", label = "Filter Data"),
                                         br(),
                                         actionButton("select_genes_bttn", label = "View Median Expression"),
                                         br(),
                                         actionButton("get_gtex_dbs", label = "Get detailed data")
                            ), 
                            mainPanel(
                              fluidRow(
                              plotlyOutput("tpm_heat")
                              ), 
                              column(width=2,
                                     br(),
                                     radioGroupButtons("gene_tpm_reads", label = "TPM or read count?", 
                                                       choices = c("TPM", "Read Count"), 
                                                       direction = "vertical")
                              ),
                              column(width=10,
                                     tags$h3("Gene Expression Distribution"),
                                     withSpinner(
                                       plotlyOutput("gene_exp"), color = '#454545'
                                     )
                              ), 
                              tags$h3("Isoform/Exon/Junction Expression"),
                              tagList(
                                column(offset = 1, width = 10,
                                       uiOutput("title")),
                                genevis_ui("gtex_plot")
                              )
                            )
                          )
                        )
                        
               )
)
               
               
               