###### Clin SeqR Shiny app ######

# this is a two file shinly application the other
# one is called server.R

# Alper Celik




server<-function(input, output, session){
  
  library(optparse)
  option_list <- list( 
    make_option(c("-c", "--credentials" , action="store", default=NULL,
                  help="Print extra output [default]")))
  
  #credentials<-unlist(strsplit(args$credentials, split = ":"))
  credentials<-c("alper", "pass")
  conn<-dbConnect(drv = RPostgreSQL::PostgreSQL(), host="localhost", 
                  dbname="gtex",
                  user=credentials[1], password=credentials[2])
  
  toastr_info(message = "loading dependencies please wait...")
  
  query<-"select geneid,genesymbol  from annotation.median_expression"
  query<-sqlInterpolate(conn, query)
  gene_table<-dbGetQuery(conn, query)
  
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
  suppressPackageStartupMessages(library(purrr))
  
  options(warn=-1)
  
  toastr_info(message = "gathering resources...")
  
  source("../modules/geneviz.R")
  source("../modules/sample_subject_filter.R")
  source("../utils/getdata.R")
  
  output$tissue_select_ui<-renderUI({
    query<-"select distinct(smtsd) from samples.samples;"
    query<-sqlInterpolate(conn, query)
    choices<-dbGetQuery(conn, query)[,1]
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
      query<-sqlInterpolate(conn, query)
      panels<-dbGetQuery(conn, query)
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
      genes<-dbGetQuery(conn, query)
    } else if (input$gtex_selection_type == "Enter Text"){
      typed<-unlist(strsplit(tm::stripWhitespace(input$gtex_text), " "))
      found<-which(gene_table$geneid %in% typed)
      not_found<-typed[!(typed %in% gene_table$geneid)]
      if(length(not_found)>0){
        toastr_warning(message = paste(not_found, "not in the gene list\n"))
      }
      gen<-paste(paste0("'", gene_table$geneid[found], "'"), collapse=",")
      query<-paste0("select * from annotation.median_expression where geneid in (", gen, ")")
      genes<-dbGetQuery(conn, query)
    } else {
      query<-"select genesymbol from annotation.gene_panels where panelname = ?panel"
      panel<-gsub(" ", "_", input$gene_panel_selection)
      query<-sqlInterpolate(conn, query, panel=panel)
      panel_genes<-unlist(dbGetQuery(conn, query))
      panel_genes<-paste("'", panel_genes, "'", sep="", collapse = ",")
      query<- paste("select * from annotation.median_expression where genesymbol in (", panel_genes, ")", sep = "")
      genes<-dbGetQuery(conn, query)
    }
    if(nrow(genes)>50 && input$gtex_selection_type!="Select Gene Panel"){
      toastr_error("Too many genes selected limit 50")
      return(NULL)
    } else if (nrow(genes)==0){
      toastr_error("No genes selected")
      return(NULL)
    } else {
      df<-list(genes=genes)
      return(df)
    }
  })
  
  output$tpm_heat<-renderPlotly({
    if(!(is.null(genes_df()))){
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
      toastr_error("Please select some genes to display their median expression values")
    }
  })
  
  tissues<-reactive({input$selected_tissues})
  
  filtered_samples<-callModule(module = filter_modal_server, id="gtex_filter", 
                               tissues=tissues, conn=conn)
  
  genes_data<-eventReactive(c(input$get_gtex_dbs,input$gene_tpm_reads, input$apply_filters), {
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
      genes_data<-get_expression(conn=conn, 
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
      toastr_warning("Please select tissues to see their expression patterns")
      NULL
    } else {
      toastr_info("click on a boxplot to see Isoform, Exon and Junction data")
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
  
  clicked_gene<-isolate({
    clicked_gene
  })
  
  callModule(module = genevis, id = "gtex_plot", 
             tissues=tissues, samples=filtered_samples,
             conn=conn, gene_id=clicked_gene)
  
  
  
  session$onSessionEnded(
    function(){
      dbDisconnect(conn)
    }
  )
}







