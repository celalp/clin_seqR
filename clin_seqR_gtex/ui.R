###### Clin SeqR Shiny app ######

# this is a two file shinly application the other
# one is called server.R

# Alper Celik



suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(shinydashboard))
suppressPackageStartupMessages(library(plotly))
suppressPackageStartupMessages(library(DT))
suppressPackageStartupMessages(library(shinytoastr))
suppressPackageStartupMessages(library(shinyWidgets))
suppressPackageStartupMessages(library(shinycssloaders))
suppressPackageStartupMessages(library(shinydashboardPlus))
suppressPackageStartupMessages(library(shinytoastr))
suppressPackageStartupMessages(library(shinyjs))
suppressPackageStartupMessages(library(DBI))
suppressPackageStartupMessages(library(RSQLite))
suppressPackageStartupMessages(library(dashboardthemes))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(shinyBS))


source("../modules/geneviz.R")
source("../modules/sample_subject_filter.R")

median_connect<-dbConnect(drv=RSQLite::SQLite(), dbname="../data/namedb.db")

ui<-dashboardPagePlus(
  dashboardHeaderPlus(enable_rightsidebar = F, title = "Clin SeqR"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Explore GTEx", tabName = "gtex", selected = T)
    )
  ),
  dashboardBody(
    fluidPage(
      shinyDashboardThemes(theme = "grey_light"),
      useToastr(),
      useShinyjs(),
      tabItems(
        tabItem(
          tabName = "gtex",
          fluidRow(
            box(title = "Select genes", width = 12, collapsible = T, collapsed = F, 
                column(width=6,
                       radioGroupButtons(inputId = "gtex_selection_type", 
                                         label = "", 
                                         choices = c("Select From List", "Enter Text"),#, "Select Gene Panel"), 
                                         justified = TRUE,
                                         selected = "Enter Text",individual = F, direction = "vertical")
                ),
                column(width=6, 
                       withSpinner(
                         uiOutput("gene_select"), color = '#8D8D8D'
                       ), 
                       br(),
                       actionButton("select_genes_bttn", label = "Select")
                )
            ),
            box(title = "Median Gene Expression", width = 12, 
                collapsible = T, collapsed = F, 
                column(width=2, 
                       pickerInput("selected_tissues", label = "Select Tissues",
                                   choices = dbGetQuery(
                                     median_connect, "PRAGMA table_info(median_expression)")[2][-c(1:2),], 
                                   options=list(`max-options`=10, 
                                                `live-search`=T, size=10), 
                                   multiple = TRUE), 
                       actionButton("get_gtex_dbs", label = "Get detailed data"), 
                       bsModal("sample_filter_modal", title = "Filter GTEx Data", 
                               trigger = "filter_gtex_data", size = "large", 
                               filter_modal_ui("gtex_filter")),
                       actionButton("filter_gtex_data", label = "Filter Data")
                ),
                column(width=10,
                       plotlyOutput("tpm_heat")
                )
            ),
            box(title="Gene Expression Distribution", width=12, 
                collapsible = T, collapsed = F, 
                column(width=2,
                       br(),
                       radioGroupButtons("gene_tpm_reads", label = "TPM or read count?", 
                                         choices = c("TPM", "Read Count"), 
                                         direction = "vertical")
                ),
                column(width=10,
                       withSpinner(
                         plotlyOutput("gene_exp"), color = '#454545'
                       )
                )
            ), 
            box(title="Isoform/Exon/Junction Expression", width = 12, 
                collapsible = T, collapsed = F, 
                tagList(
                  column(offset = 1, width = 10,
                         uiOutput("title")),
                  genevis_ui("gtex_plot")
                )
            )
          )
        )
      )
      #,
      #tabItem(
      #  tabName = "expression"
      # put boxes
      #), 
      #tabItem(
      # tabName = "variants"
      #variant module here
      #)
    )
  )
)


