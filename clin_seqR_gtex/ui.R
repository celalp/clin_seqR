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
suppressPackageStartupMessages(library(shinyjs))
suppressPackageStartupMessages(library(DBI))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(shinyBS))
suppressPackageStartupMessages(library(shinythemes))
suppressPackageStartupMessages(library(shinydashboard))
suppressPackageStartupMessages(library(dashboardthemes))


source("../modules/geneviz.R")
source("../modules/sample_subject_filter.R")
source("../utils/getdata.R")

ui<-dashboardPage(
  dashboardHeader(title = "ClinSeqR"), 
  dashboardSidebar(#put icons
    sidebarMenu(
      menuItem("Explore GTEx Data", tabName="gtex", icon=icon("dashboard")), 
      menuItem("Samples", icon = icon("database"), startExpanded = F,
               menuSubItem("Expression", tabName = "expression", icon=icon("bar-chart")),
               menuSubItem("Variants", tabName = "variants", icon=icon("table"))
      ), 
      menuItem("Admin Console", tabName="admin", icon=icon("search"))
    ),
    actionButton("login_button", label = "Login", icon=icon("user"))
  ), 
  dashboardBody(
    shinyDashboardThemes(theme = "grey_light"),
    useToastr(),
    tabItems(
      tabItem(tabName = "gtex",
        fluidRow(
          box(title = "Select Gene(s)", width = 4,
              radioGroupButtons(inputId = "gtex_selection_type", 
                                label = "", 
                                choices = c("Select From List", "Enter Text", "Select Gene Panel"), 
                                justified = TRUE, status = "primary",
                                selected = "Select From List",individual = F, direction = "vertical"),
              uiOutput("gene_select"), color = '#8D8D8D',
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
          box(title = "Median Expression", width = 8, 
                tagList(
                  bsAlert("genes_alert_heat"),
                  plotlyOutput("tpm_heat", height = "600px")
                )
              )
        ),
        fluidRow(
          box(title = "Gene Expression Distribution", width = 12,
              column(width=2,
                     br(),
                     radioGroupButtons("gene_tpm_reads", label = "TPM or read count?", 
                                       choices = c("TPM", "Read Count"), 
                                       direction = "vertical")
              ),
              column(width=10,
                     tagList(
                       bsAlert("genes_alert_box"),
                       withSpinner(
                         plotlyOutput("gene_exp"), color = '#454545'
                       ))
              )
          ),
          box(title = "Isoform/Exon/Junction Expression",  width = 12, collapsible = T, 
              collapsed = T,
              tagList(
                uiOutput("title"),
                genevis_ui("gtex_plot")
              )
          ))
      ), 
      tabItem(tabName = "expression"), 
      
      tabItem(tabName = "variants"), 
      tabItem(tabName = "admin")
    )
  )
)
    
    
    
    