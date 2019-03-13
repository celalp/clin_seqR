
gene_select_ui<-function(id){
  ns<-NS(id)
  tagList(
    uiOutput(ns("selection_mode")),
    uiOutput(ns("gene_select")),
    br(),
    uiOutput(ns("buttons"))
  )
}


gene_select_server<-function(input, output, session, conn, login){
  output$buttons<-renderUI({
    if(login$login){
      tagList(
        actionButton(session$ns("get_conn_dbs"), label = "Select Genes")
      )
    } else {
      NULL
    }
  })
  
  output$selection_mode<-renderUI({
    if(login$login){
      radioGroupButtons(inputId = session$ns("conn_selection_type"), 
                        label = "", 
                        choices = c("Select From List", "Enter Text", "Select Gene Panel"), 
                        justified = TRUE, status = "primary",
                        selected = "Select From List",individual = F, direction = "vertical")
    } else {
      NULL
    }
  })
  
  query<-"select geneid,genesymbol from annotation.gene_names"
  query<-sqlInterpolate(conn, query)
  gene_table<-dbGetQuery(conn, query)
  
  output$gene_select<-renderUI({
    if(login$login){
      if(input$conn_selection_type=="Select From List"){
        output$selection_table<-renderDT({
          datatable(gene_table, selection = 'multiple', class = 'compact', 
                    options=list(Dom='f', extensions = c('Responsive'), paging=T, 
                                 deferRender=T), rownames = F)}, server = T)
        tagList(
          DTOutput(session$ns("selection_table")),
          actionButton(session$ns("deselect"), label = "Deselect")
        )
      } else if (input$conn_selection_type=="Enter Text"){
        textAreaInput(session$ns("conn_text"), label = "Enter genes", resize = 'vertical',
                      placeholder = "Enter Ensembl IDs here one per line")
      } else {
        query<-"select distinct(gene_panels.panelname) from annotation.gene_panels"
        query<-sqlInterpolate(conn, query)
        panels<-dbGetQuery(conn, query)
        panels<-gsub("_", " ", panels$panelname)
        pickerInput(session$ns("gene_panel_selection"), "Select Gene Panel", multiple = F, 
                    choices = panels)
      }
    } else {
      NULL
    }
  })
  
  gene_select_proxy<-dataTableProxy("selection_table", session = session)
  observeEvent(input$deselect, {
    selectRows(gene_select_proxy, NULL)
  })
  
  genes<-eventReactive(input$get_conn_dbs, {
    if(input$conn_selection_type=="Select From List"){
      genes<-gene_table[input$selection_table_rows_selected,]
    } else if (input$conn_selection_type == "Enter Text"){
      typed<-unlist(strsplit(tm::stripWhitespace(input$conn_text), " "))
      found<-which(gene_table$geneid %in% typed)
      not_found<-typed[!(typed %in% gene_table$geneid)]
      if(length(not_found)>0){
        toastr_warning(message = paste(not_found, "not in the gene list\n"))
      }
      genes<-gene_table[found,]
    } else {
      query<-"select genesymbol from annotation.gene_panels where panelname = ?panel"
      panel<-gsub(" ", "_", input$gene_panel_selection)
      query<-sqlInterpolate(conn, query, panel=panel)
      panel_genes<-unlist(dbGetQuery(conn, query))
      genes<-gene_table %>%
        filter(genesymbol %in% panel_genes)
    }
    return(genes)
  })
  
  return(genes)
}