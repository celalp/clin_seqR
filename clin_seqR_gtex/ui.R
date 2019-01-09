###### Clin SeqR Shiny app ######

# this is a two file shinly application the other
# one is called server.R

# Alper Celik

#TODO add bsalerts to all tabs including the isoform one and remove corresponding toastr alerts
#TODO same for the filter data ui
#TODO collapsible gene expression and isoform module boxes to be expanded with some thing (actionlink?)

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
source("../utils/getdata.R")

ui<-navbarPageWithInputs("GTEx Module", inverse = F, theme = shinytheme("cerulean"), 
                         inputs = actionButton("login_modal", "login", style="margin:0.5%"),
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
                                                 tagList(
                                                   filter_modal_ui("gtex_filter"),
                                                   actionButton("apply_filters", "Apply Filters")
                                                 )),
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
                                  bsAlert("genes_alert_heat"),
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
                                     bsAlert("genes_alert_box"),
                                     withSpinner(
                                       plotlyOutput("gene_exp"), color = '#454545'
                                     )
                              ), 
                              tags$h3("Isoform/Exon/Junction Expression"),
                              tagList(
                                uiOutput("title"),
                                genevis_ui("gtex_plot")
                              )
                            )
                          )
                        )
               ),
               tabPanel("Sample Expression"),
               tabPanel("Sample Variants"), 
               tabPanel("Admin Console")
)


