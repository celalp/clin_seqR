###### Clin SeqR Shiny app ######

# this is a two file shinly application the other
# one is called server.R

# Alper Celik




server<-function(input, output, session){
  
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
  suppressPackageStartupMessages(library(shinyjs))
  suppressPackageStartupMessages(library(DBI))
  suppressPackageStartupMessages(library(reshape2))
  suppressPackageStartupMessages(library(viridis))
  suppressPackageStartupMessages(library(tm))
  suppressPackageStartupMessages(library(data.table))
  suppressPackageStartupMessages(library(shinyBS))
  
  options(warn=-1)
  
  toastr_info(message = "gathering resources...")
  
  login<-F
  
  login_modal<-function(login){
    modalDialog(title = "Login to database",
    textInput(inputId = "username", label = "Username"),
    passwordInput(inputId = "password", label = "Password"),
    actionButton(inputId = "modal_login", label = "Login", icon=icon("user"))
    )
  }
  
  observeEvent(input$login_button, {
    showModal(login_modal())
  })
  
  #observeEvent()
  
  source("../modules/geneviz.R")
  source("../modules/sample_subject_filter.R")
  source("../utils/getdata.R")
  
  #login<-eventReactive()
  # will create a new user with select priviliges to specific schemas to login here
  # other tabs will be username/password specific
  #credentials<-c("gtexuser", "gtexuserforclinseqr")
  credentials<-c("alper", "pass")
  gtex<-dbConnect(drv = RPostgreSQL::PostgreSQL(), host="localhost", 
                  dbname="gtex",
                  user=credentials[1], password=credentials[2])
  
  query<-"select geneid,genesymbol from annotation.median_expression"
  query<-sqlInterpolate(gtex, query)
  gene_table<-dbGetQuery(gtex, query)
  
  
  #####################################################
  ########### gtex tab server #########################
  #####################################################
  
  output$tissue_select_ui<-renderUI({
    query<-"select distinct(smtsd) from samples.samples;"
    query<-sqlInterpolate(gtex, query)
    choices<-dbGetQuery(gtex, query)[,1]
    choices<-unlist(gsub("_", " ", choices))
    choices<-choices[order(choices)]
    pickerInput("selected_tissues", label = "Select Tissues",
                choices = choices, 
                options=list(`max-options`=10, 
                             `live-search`=T, size=10), 
                multiple = TRUE)
  })
  
  output$gene_select<-renderUI({
    if(input$gtex_selection_type=="Select From List"){
      output$selection_table<-renderDT({
        datatable(gene_table, selection = 'multiple', class = 'compact', 
                  options=list(Dom='f', extensions = c('Responsive'), paging=T, 
                               deferRender=T), rownames = F)}, server = T)
      tagList(
        DTOutput("selection_table"),
        actionButton("deselect", label = "Deselect")
      )
    } else if (input$gtex_selection_type=="Enter Text"){
      textAreaInput("gtex_text", label = "Enter genes", resize = 'vertical',
                    placeholder = "Enter Ensembl IDs here one per line. Max 20")
    } else {
      query<-"select distinct(gene_panels.panelname) from annotation.gene_panels"
      query<-sqlInterpolate(gtex, query)
      panels<-dbGetQuery(gtex, query)
      panels<-gsub("_", " ", panels$panelname)
      pickerInput("gene_panel_selection", "Select Gene Panel", multiple = F, 
                  choices = panels)
    }
  })
  
  gene_select_proxy<-dataTableProxy("selection_table", session = session)
  observeEvent(input$deselect, {
    selectRows(gene_select_proxy, NULL)
  })
  
  genes_df<-eventReactive(c(input$select_genes_bttn, input$get_gtex_dbs), {
    if(input$gtex_selection_type=="Select From List"){
      genes<-gene_table$geneid[input$selection_table_rows_selected]
      genes<-paste(paste0("'", genes, "'"), collapse=",")
      query<-paste("select * from annotation.median_expression where geneid in (", genes, ")")
      genes<-dbGetQuery(gtex, query)
    } else if (input$gtex_selection_type == "Enter Text"){
      typed<-unlist(strsplit(tm::stripWhitespace(input$gtex_text), " "))
      found<-which(gene_table$geneid %in% typed)
      not_found<-typed[!(typed %in% gene_table$geneid)]
      if(length(not_found)>0){
        toastr_warning(message = paste(not_found, "not in the gene list\n"))
      }
      gen<-paste(paste0("'", gene_table$geneid[found], "'"), collapse=",")
      query<-paste0("select * from annotation.median_expression where geneid in (", gen, ")")
      genes<-dbGetQuery(gtex, query)
    } else {
      query<-"select genesymbol from annotation.gene_panels where panelname = ?panel"
      panel<-gsub(" ", "_", input$gene_panel_selection)
      query<-sqlInterpolate(gtex, query, panel=panel)
      panel_genes<-unlist(dbGetQuery(gtex, query))
      panel_genes<-paste("'", panel_genes, "'", sep="", collapse = ",")
      query<- paste("select * from annotation.median_expression where genesymbol in (", panel_genes, ")", sep = "")
      genes<-dbGetQuery(gtex, query)
    }
    if(nrow(genes)>50 && input$gtex_selection_type!="Select Gene Panel"){
      createAlert(session, "genes_alert_heat", "genes_alert_heat_control", title = "",
                  content = "Too many genes selected limit 50",style = "warning", 
                  append = FALSE)
      return(NULL)
    } else if (nrow(genes)==0){
      return(NULL)
    } else {
      df<-list(genes=genes)
      return(df)
    }
  })
  
  output$tpm_heat<-renderPlotly({
    if(!(is.null(genes_df()))){
      closeAlert(session, "genes_alert_heat_control")
      forheat<-genes_df()$genes
      tissues<-colnames(forheat)[-c(1:3)]
      genes<-forheat$genesymbol
      mat<-as.matrix(forheat[,-c(1:3)])
      ## rescale colors
      vals <- unique(scales::rescale(c(mat)))
      o <- order(vals, decreasing = FALSE)
      cols <- scales::col_numeric("Spectral", domain = NULL)(vals)
      colz <- setNames(data.frame(vals[o], cols[o]), NULL)
      ####
      plot_ly(y=tissues, x=genes, z=t(mat), type="heatmap", 
              source="tpm_heatplot", colorscale=colz) %>% 
        layout(xaxis=list(title="", dtick=1), 
               yaxis=list(title="", dtick=1))
    } else {
      createAlert(session, "genes_alert_heat", "genes_alert_heat_control", title = "",
                  content = "No genes selected please select genes or gene panels from the panel on the left",style = "warning", 
                  append = FALSE)
    }
  })
  
  tissues<-reactive({input$selected_tissues})
  
  filtered_samples<-callModule(module = filter_modal_server, id="gtex_filter", 
                               tissues=tissues, conn=gtex)
  
  genes_data<-reactive({
    if(length(input$selected_tissues)==0){
      genes_data<-NULL
    } else {
      if(input$gene_tpm_reads=="TPM"){
        table="gene_tpm"
      } else {
        table="gene_reads"
      }
      gene_names<-genes_df()$genes
      tissues<-input$selected_tissues
      genes_data<-get_expression(conn=gtex, 
                                 gene_id=genes_df()$genes$geneid, 
                                 tissue=input$selected_tissues, 
                                 samples=filtered_samples(),
                                 table=table, extra_columns="gene_name"
      )
      genes_data<-do.call("rbind", genes_data)
      return(list(genes_data=genes_data))
    }
  })
  
  output$gene_exp<-renderPlotly({
    if(is.null(genes_data())){
      createAlert(session, "genes_alert_box", "genes_alert_box_control", title = "",
                  content = "No genes selected please select genes or gene panels from the panel on the left",style = "warning", 
                  append = FALSE)
      NULL
    } else {
      closeAlert(session, "genes_alert_box_control")
      createAlert(session, "genes_alert_box", "genes_alert_box_control", title = "",
                  content = "Click on the boxplot to see detailed gene expression",
                  style = "info", append = FALSE)
      plot_ly(data = genes_data()$genes_data, x=~tissue, y=~value, 
              color=~gene_name, type = "box", source = "gene_boxplot")%>%
        layout(boxmode = "group", xaxis=list(title="Tissue"), 
               yaxis=list(title=input$gene_tpm_reads)
        )
    }
  })
  
  clicked_gene<-reactive({
    s <- event_data("plotly_click", source = "gene_boxplot")
    if (length(s) == 0) {
      gene_id<-NULL
      gene_symbol<-NULL
    } else {
      # plotly sorts by gene name but returns a 0 based index so I need to sort gene names
      # then get +1 index
      if(input$gene_tpm_reads=="TPM"){
        dat_tbl<-"transcript_tpm"
      } else {
        dat_tbl<-"transcript_reads"
      }
      s<-as.list(s)
      all_gene_names<-sort(unique(genes_df()$genes$genesymbol))
      selected_gene<-unique(all_gene_names[(s$curveNumber+1)])
      gene_symbol<-genes_df()$genes$genesymbol[genes_df()$genes$genesymbol==selected_gene]
      gene_id<-genes_df()$genes$geneid[genes_df()$genes$genesymbol==selected_gene]
    }
    return(c(gene_id, gene_symbol))
  })
  
  output$title<-renderUI({
    if(!is.null(clicked_gene())){
      tags$h4(paste("Displaying information for ", clicked_gene()[2]))
    } else {
      NULL
    }
  })
  callModule(module = genevis, id = "gtex_plot", 
             tissues=tissues, samples=filtered_samples,
             conn=gtex, gene_id=clicked_gene)
  
  
  #####################################################
  ########### sample expression #######################
  #####################################################
  
  #callModule(module=sample_expression, id="sample_expression", login=login)
  
  #####################################################
  ########### sample variants #########################
  #####################################################
  
  
  #####################################################
  ############# admin console #########################
  #####################################################
  
  
  
  session$onSessionEnded(
    function(){
      dbDisconnect(gtex)
    }
  )  
}







