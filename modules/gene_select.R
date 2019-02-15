
gene_select_ui<-function(id){
  ns<-NS(id)
  tagList(
    uiOutput(ns("selection_mode")),
    uiOutput(ns("gene_select")),
    br(),
    uiOutput(ns("tissue_select_ui")),
    uiOutput(ns("buttons"))
  )
}


gene_select<-function(input, output, session, conn, login){
  output$buttons<-renderUI({
    if(login$login){
      tagList(
        actionButton(session$ns("filter_conn_data"), label = "Filter Data"),
        br(),
        actionButton(session$ns("get_conn_dbs"), label = "Get Detailed Data")
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
  
  output$tissue_select_ui<-renderUI({
    if(login$login){
      query<-"select distinct(smtsd) from samples.samples;"
      query<-sqlInterpolate(conn, query)
      choices<-dbGetQuery(conn, query)[,1]
      choices<-unlist(gsub("_", " ", choices))
      choices<-choices[order(choices)]
      pickerInput(session$ns("selected_tissues"), label = "Select Tissues",
                  choices = choices, 
                  options=list(`max-options`=10, 
                               `live-search`=T, size=10), 
                  multiple = TRUE)
    } else {
      NULL
    }
  })
  
  output$gene_select<-renderUI({
    if(login$login){
      query<-"select geneid,genesymbol from annotation.median_expression"
      query<-sqlInterpolate(conn, query)
      gene_table<-dbGetQuery(conn, query)
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
                      placeholder = "Enter Ensembl IDs here one per line. Max 20")
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
  
  observeEvent(input$filter_conn_data, {
    showModal(
      modalDialog(size="l", easyClose = T, footer = modalButton("Close"),
                  tagList(
                  filter_modal_ui(session$ns("gtex_filter")),
                  actionButton(session$ns("apply_filters"), "Apply Filters")
                  )
    ))
  })
  
  tissues<-reactive({input$selected_tissues})
  
  filtered_samples<-callModule(module = filter_modal_server, id="gtex_filter", 
             tissues=tissues, conn=conn)
  
  genes_df<-eventReactive(input$get_gtex_dbs, {
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
  
  #TODO return a proper list here for the other module to handle
  return(
    reactive({list(tissues=tissues())})
    )#, samples=filtered_samples)), genes_df=isolate(genes_df())))
  
}