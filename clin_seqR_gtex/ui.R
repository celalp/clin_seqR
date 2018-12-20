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
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(shinyBS))
suppressPackageStartupMessages(library(shinythemes))


source("../modules/geneviz.R")
source("../modules/sample_subject_filter.R")

ui<-navbarPage("GTEx Module", inverse = T, theme = shinytheme("cerulean"),
               tabPanel("Explore GTEx Data", 
                        fluidRow(
                          sidebarLayout(
                            sidebarPanel(width=3,
                                         tags$h3("Select Gene(s)"),
                                         radioGroupButtons(inputId = "gtex_selection_type", 
                                                           label = "", 
                                                           choices = c("Select From List", "Enter Text", "Select Gene Panel"), 
                                                           justified = TRUE, status = "primary",
                                                           selected = "Select From List",individual = F, direction = "vertical"),
                                         withSpinner(
                                           uiOutput("gene_select"), color = '#8D8D8D'
                                         ), 
                                         br(),
                                         uiOutput("tissue_select_ui"),
                                         bsModal("sample_filter_modal", title = "Filter GTEx Data", 
                                                 trigger = "filter_gtex_data", size = "large", 
                                                 filter_modal_ui("gtex_filter")),
                                         actionButton("select_genes_bttn", label = "Display Median Expression"),
                                         br(),
                                         actionButton("filter_gtex_data", label = "Filter Data"),
                                         br(),
                                         actionButton("get_gtex_dbs", label = "Get Detailed Data")
                            ), 
                            mainPanel(
                              useToastr(),
                              fluidRow(
                                tagList(
                                  #this does not display properly
                                tags$h3("Median Expression"),
                                plotlyOutput("tpm_heat", height = "600px")
                                )
                              ), 
                              column(width=2,
                                     br(),
                                     radioGroupButtons("gene_tpm_reads", label = "TPM or read count?", 
                                                       choices = c("TPM", "Read Count"), 
                                                       direction = "vertical")
                              ),
                              column(width=10,
                                     tags$h3("Gene Expression Distribution"),
                                     #bsAlert("gene_exp_alert"),
                                     withSpinner(
                                       plotlyOutput("gene_exp"), color = '#454545'
                                     )
                              ), 
                              tags$h3("Isoform/Exon/Junction Expression"),
                              tagList(
                                #bsAlert("genevis_alert"),
                                uiOutput("title"),
                                genevis_ui("gtex_plot")
                              )
                            )
                          )
                        )
                        
               )
)


