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
  suppressPackageStartupMessages(library(shinytoastr))
  suppressPackageStartupMessages(library(shinyjs))
  suppressPackageStartupMessages(library(DBI))
  suppressPackageStartupMessages(library(reshape2))
  suppressPackageStartupMessages(library(GenomicFeatures))
  suppressPackageStartupMessages(library(viridis))
  suppressPackageStartupMessages(library(tm))
  suppressPackageStartupMessages(library(data.table))
  suppressPackageStartupMessages(library(RSQLite))
  suppressPackageStartupMessages(library(shinyBS))
  suppressPackageStartupMessages(library(purrr))
  
  options(warn=-1)
  
  toastr_info(message = "gathering resources...")
  
  datadir<-"../data/"
  source("../modules/geneviz.R")
  source("../modules/sample_subject_filter.R")
  source("../utils/getdata.R")
  median_connect<-dbConnect(drv=RSQLite::SQLite(), dbname="../data/namedb.db")
  gene_table<-dbGetQuery(median_connect,"select [Ensembl Id], [Gene Name] from median_expression")
  annots<-loadDb("../data/gtex_annotation.db")
  collapsed_annots<-loadDb("../data/gtex_collapsed_annotation.db")
  removeClass(selector = "body", class = "sidebar-collapse")
  ######## GTEX TAB I/O #########
  
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
      textAreaInput("gtex_text", label = "Enter genes", 
                    placeholder = "Enter Ensembl IDs here one per line. Max 20")
    } else {
      panels<-dbGetQuery(median_connect, "select distinct(panel_name) from gene_panels")
      pickerInput("gene_panel_selection", "Select Gene Panel", multiple = F, 
                  choices = panels)
    }
  })
  
  gene_select_proxy<-dataTableProxy("selection_table", session = session)
  observeEvent(input$deselect, {
    selectRows(gene_select_proxy, NULL)
  })
  
  genes_df<-eventReactive(input$select_genes_bttn, {
    if(input$gtex_selection_type=="Select From List"){
      rows<-paste(input$selection_table_rows_selected, collapse = ",")
      genes<-dbGetQuery(median_connect, 
                        paste("select * from median_expression where rowid in (", rows, ")", collapse = ""))
    } else if (input$gtex_selection_type == "Enter Text"){
      typed<-unlist(strsplit(input$gtex_text, "\n"))
      found<-which(gene_table[,1] %in% typed)
      genes<-dbGetQuery(median_connect, 
                        paste("select * from median_expression where rowid in (", rows, ")", collapse = ""))
      not_found<-typed[!(typed %in% genes[,1])]
      if(length(not_found)>0){
        toastr_warning(message = paste(not_found, "not in the gene list\n"))
      }
    } else {
      panel_genes<-unlist(dbGetQuery(median_connect, paste("select gene_symbol from gene_panels where panel_name = '", 
                                                           input$gene_panel_selection, "'", sep = "")))
      panel_genes<-paste("'", panel_genes, "'", sep="", collapse = ",")
      query<- paste("select * from median_expression where [Gene Name] in (", panel_genes, ")", sep = "")
      genes<-dbGetQuery(median_connect, query)
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
      tissues<-colnames(forheat)[-c(1:2)]
      genes<-forheat$'Gene Name'
      mat<-as.matrix(forheat[,-c(1:2)])
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
    }
  })
  
  tissues<-reactive({input$selected_tissues})
  
  filtered_samples<-callModule(module = filter_modal_server, id="gtex_filter", 
                               tissues=tissues, datadir=datadir)
  
  genes_data<-eventReactive(c(input$get_gtex_dbs,input$gene_tpm_reads), {
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
      genes_data<-get_expression(datadir=datadir, 
                                 gene_id=genes_df()$genes$'Ensembl Id', 
                                 tissue=input$selected_tissues, 
                                 samples=filtered_samples(), # this will be samples()$samples
                                 table=table, extra_columns="gene_name"
      )
      genes_data<-do.call("rbind", genes_data)
      return(list(genes_data=genes_data))
    }
  })
  
  output$gene_exp<-renderPlotly({
    if(is.null(genes_data()$genes_data) && !is.null(genes_df())){
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
      expression_data<-NULL
    } else {
      # plotly sorts by gene name but returns a 0 based index so I need to sort gene names
      # then get +1 index
      if(input$gene_tpm_reads=="TPM"){
        dat_tbl<-"transcript_tpm"
      } else {
        dat_tbl<-"transcript_reads"
      }
      s<-as.list(s)
      all_gene_names<-sort(unique(genes_df()$genes$`Gene Name`))
      selected_gene<-unique(all_gene_names[(s$curveNumber+1)])
      gene_id<-genes_df()$genes$`Ensembl Id`[genes_df()$genes$`Gene Name`==selected_gene]
    }
    return(gene_id)
  })
  
  output$title<-renderUI({
    if(!is.null(clicked_gene())){
        tags$h4(paste("Displaying information for ", clicked_gene()))
    } else {
      NULL
    }
  })
  
  clicked_gene<-isolate({
    clicked_gene
  })
  
  callModule(module = genevis, id = "gtex_plot", 
             datadir=datadir, tissues=tissues, samples=filtered_samples,
             annot=annots, collapsed_annot=collapsed_annots,
             gene_id=clicked_gene)

}







