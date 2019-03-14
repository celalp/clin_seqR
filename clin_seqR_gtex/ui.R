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
suppressPackageStartupMessages(library(shinydashboard))
suppressPackageStartupMessages(library(dashboardthemes))
suppressPackageStartupMessages(library(shinydashboardPlus))


#TODO include conditional panels for alerting when selections are not valid

source("../modules/geneviz.R")
source("../utils/getdata.R")
source("../modules/gene_select.R")
source("../modules/tissue_select.R")
source("../modules/gene_expression.R")

ui<-dashboardPagePlus(
  dashboardHeaderPlus(title = "ClinSeqR", 
                      left_menu = tagList(
                        dropdownBlock(
                          id = "loginmenu",
                          icon = "user", title = "Login/Logout",
                          uiOutput("login_menu"), badgeStatus = NULL
                        )
                      )
  ),
  dashboardSidebar(#put icons
    sidebarMenu(
      menuItem("Explore GTEx Data", tabName="gtex", icon=icon("dashboard")), 
      menuItem("Samples", icon = icon("database"), startExpanded = F,
               menuSubItem("Expression", tabName = "expression", icon=icon("bar-chart")),
               menuSubItem("Variants", tabName = "variants", icon=icon("table"))
      ), 
      menuItem("Account Settings", tabName="account", icon=icon("user"))
    )
  ), 
  dashboardBody(
    fluidRow(
      shinyDashboardThemes(theme = "grey_light"),
      useToastr(), 
      useShinyjs(),
      uiOutput("login_prompt"),
      tabItems(
        tabItem(tabName = "gtex", #TODO need to add bs alerts to median expression, gene expression dist, and genevis
                box(title = "Select Gene(s)", width = 4, 
                    gene_select_ui("gene_select")), 
                box(title = "Select tissue(s)", width = 8,
                    tissue_select_ui("tissue_select")),
                box(title = "Gene Expression Distribution", width = 12, collapsible = T,
                    gene_expression_ui("gene_expression")),
                box(title = "Isoform/Exon/Junction Expression",  width = 12, collapsible = T, 
                    collapsed = F,
                    genevis_ui("detailed_expression"))
        ), 
        tabItem(tabName = "expression"), 
        tabItem(tabName = "variants")
      )
    )
  )
)
