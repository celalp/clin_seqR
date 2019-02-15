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
suppressPackageStartupMessages(library(shinydashboard))
suppressPackageStartupMessages(library(dashboardthemes))
suppressPackageStartupMessages(library(shinydashboardPlus))


source("../modules/geneviz.R")
source("../modules/sample_subject_filter.R")
source("../utils/getdata.R")
source("../modules/gene_select.R")
source("../modules/median_heatmap.R")

ui<-dashboardPagePlus(
  dashboardHeaderPlus(title = "ClinSeqR", 
                      left_menu = tagList(
                        dropdownBlock(
                          id = "loginmeni",
                          title = "Login",
                          icon = icon("user"),
                          tagList(
                            textInput(inputId = "username", label = "Username"),
                            passwordInput(inputId = "password", label = "Password"),
                            actionButton(inputId = "login_button", label = "Login", icon=icon("user"))
                          ), badgeStatus = NULL
                        )
                      )
  ),
  dashboardSidebar(#put icons
    sidebarMenu(
      menuItem("Explore GTEx Data", tabName="gtex", icon=icon("dashboard")), 
      menuItem("Samples", icon = icon("database"), startExpanded = F,
               menuSubItem("Expression", tabName = "expression", icon=icon("bar-chart")),
               menuSubItem("Variants", tabName = "variants", icon=icon("table"))
      ) 
      #menuItem("Admin Console", tabName="admin", icon=icon("search"))
    )
  ), 
  dashboardBody(
    shinyDashboardThemes(theme = "grey_light"),
    useToastr(), 
    useShinyjs(), 
    tabItems(
      tabItem(tabName = "gtex",
              fluidRow(
                box(title = "Select Gene(s)", width = 4, 
                    gene_select_ui("gene_select")), 
                box(title = "Median Expression", width = 8, 
                    verbatimTextOutput("tissue_select_check"))
                #median_heatmap_ui("median_heatmap"))
              ),
              fluidRow(
                box(title = "Gene Expression Distribution", width = 12),
                box(title = "Isoform/Exon/Junction Expression",  width = 12, collapsible = T, 
                    collapsed = F,
                    tagList(
                      #uiOutput("title"),
                      genevis_ui("gtex_plot")
                    )
                )
              )
      )
    ), 
    tabItem(tabName = "expression"), 
    tabItem(tabName = "variants")
  )
  
)
