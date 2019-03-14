
genevis_ui<-function(id){
  ns<-NS(id)
  tagList(
    uiOutput(ns("gene_picker_ui")),
    tags$br(),
    tabsetPanel( 
      tabPanel("Isoform Expression", value = "gtex_isoform",
               tagList(
                 column(width=2, 
                        uiOutput(ns("selection_ui"))),
                 column(width=10,
                        tagList(
                          tags$br(),
                          uiOutput(ns("isoform_alert")),
                          plotlyOutput(ns("isoform_exp")),
                          plotOutput(ns("isoform_model"))
                        )
                 ))
      ),
      tabPanel("Exon Expression", value = "gtex_exon", 
               tagList(
                 column(width=4,
                        tags$br(),
                        dataTableOutput(ns("exon_table")),
                        uiOutput(ns("deselect_exon")),
                        tags$br(),
                        uiOutput(ns("exon_alert"))
                 ),
                 column(width=8,
                        tags$br(),
                        plotlyOutput(ns("exon_exp")),
                        tags$br(),
                        plotOutput(ns("gene_model_collapsed")))
               )
      ),
      tabPanel("Junction Expression", value="gtex_junc", 
               tagList(
                 column(width=4, 
                        tags$br(), 
                        dataTableOutput(ns("junc_table")),
                        uiOutput(ns("deselect_junc")),
                        tags$br(),
                        uiOutput(ns("junc_alert"))
                 ),
                 column(width=8, 
                        tags$br(),
                        plotlyOutput(ns("junc_exp")),
                        tags$br(),
                        plotOutput(ns("junc_gene_model"))
                 )
               ))
    )
  )
}

genevis<-function(input, output, session, tissues_samples, genes, login, conn){ # this needs data as well
  
  ready<-reactive({
    if(is.null(genes()) | is.null(tissues_samples()) | !(login$login)){
      ready<-F
    } else {
      ready<-T
    }
    return(ready)
  })
  
  txdf<-reactive({
    if(ready()){
      gene_name<-genes()$geneid[genes()$genesymbol==input$gene_picker]
      txdf_full<-parse_exon_df(conn = conn, gene_name = gene_name, collapsed=F)
      txdf_collapsed<-parse_exon_df(conn = conn, gene_name = gene_name, collapsed=T)
      txdf<-list(gene_id=gene_name, full=txdf_full, collapsed=txdf_collapsed)
    } else {
      txdf<-NULL
    }
    return(txdf)
  })
  
  expression_data<-reactive({
    if(ready()){
      if(input$isoform_tpm_reads=="TPM"){
        iso_tbl<-"transcript_tpm"
      } else {
        iso_tbl<-"transcript_reads"
      }
      isoform_expression<-get_expression(conn=conn, gene_id=txdf()$gene_id, 
                                         tissues_samples=tissues_samples(), 
                                         table=iso_tbl, extra_columns='"transcript_id"')
      exon_expression<-get_expression(conn=conn, gene_id=txdf()$gene_id, 
                                      tissues_samples=tissues_samples(), 
                                      table="exon_reads", extra_columns='"exon_id"')
      junction_expression<-get_expression(conn=conn, gene_id=txdf()$gene_id, 
                                          tissues_samples=tissues_samples(), 
                                          table="junction_reads", extra_columns=c('"junction_id"', '"start"','"end"'))
      expression_data<-list(isoform=isoform_expression, exon=exon_expression, junction=junction_expression)
    } else {
      expression_data<-NULL
    }
    return(expression_data)
  })
  
  
  ##### Isoform related #########
  
  output$isoform_alert<-renderUI({
    alert_box(alert = "Please select genes and tissues above to see expression distribution", 
              color = "maroon", condition = !(ready()))
  })
  
  output$gene_picker_ui<-renderUI({
    if(ready()){
      tagList(
        pickerInput(inputId = session$ns("gene_picker"), label = "Select a gene to view expression details", 
                    choices = genes()$genesymbol, multiple = F, options = list(`live-search`=T, size=10), 
                    selected = NULL)
      )
    } else {
      NULL
    }
  })
  
  output$selection_ui<-renderUI({
    if(ready()){
      tagList(
        tags$br(),
        radioGroupButtons(session$ns("isoform_tpm_reads"), label = "TPM or read count?", 
                          choices = c("TPM", "Read Count"), direction = "vertical", 
                          selected = "TPM"))
    } else {
      NULL
    }
  })
  
  output$isoform_exp<-renderPlotly({
    if(ready()){
      plot_ly(data = expression_data()$isoform, x=~tissue, y=~value, 
              color=~transcript_id, type = "box", source = "isoform_boxplot")%>%
        layout(boxmode = "group", xaxis=list(title="Tissue"), 
               yaxis=list(title=input$isoform_tpm_reads))
    } else {
      NULL
    }
  })
  
  output$isoform_model<-renderPlot({
    if(ready()){
      hover<-event_data(event = "plotly_hover", source="isoform_boxplot")$curveNumber+1
      p<-plot_transcripts(tx_df= txdf()$full, hover=hover)
      p
    } else {
      NULL
    }
  },height = function() {
    if(ready()){
      max(txdf()$full$tx_rank*20)
    } else {
      1
    }
  })
  
  ########### exon related ########
  
  output$exon_table<-renderDataTable({
    if(ready()){
      dat<-txdf()$collapsed
      dat<-dat[dat$class=="exon",]
      dat<-dat[order(as.integer(dat$exonrank)),]
      rownames(dat)<-c(1:nrow(dat))
      datatable(dat[,c("exonname", "exonrank")])
    } else {
      NULL
    }
  })
  
  output$deselect_exon<-renderUI({
    if(ready()){
      actionButton(session$ns("deselect_exon_bttn"), label = "Deselect")
    } else {
      NULL
    }
  })
  
  exon_select_proxy<-dataTableProxy("exon_table", session = session)
  observeEvent(input$deselect_exon_bttn, {
    selectRows(exon_select_proxy, NULL)
  })
  
  selected_exons<-reactive({
    nointron<-txdf()$collapsed[txdf()$collapsed$class=="exon",]
    exons<-nointron$exonname[input$exon_table_rows_selected]
    return(exons)
  })
  
  output$exon_alert<-renderUI({
    if(!ready()){
      msg<-"Please select genes and tissues above to see expression distribution"
      color<-"maroon" 
      condition<-T
    } else if (ready() && length(selected_exons())<1){
      msg<-"Select exons from the table above to see their expression"
      color<-"light-blue"
      condition<-T
    } else {
      msg<-""
      color<-"red"
      condition<-F
    }
    tagList(
      alert_box(alert = msg, color = color, condition = condition)
    )
  })
  
  output$exon_exp<-renderPlotly({
    if(ready()){
      plot_data<-expression_data()$exon[(expression_data()$exon$exon_id %in% selected_exons()), ]
      plot_ly(data = plot_data,
              x=~tissue, y=~value, 
              color=~exon_id, type = "box", source = "exon_boxplot")%>%
        layout(boxmode = "group", xaxis=list(title="Tissue"), 
               yaxis=list(title="Read Count"))
    }else {
      NULL
    }
  })
  
  output$gene_model_collapsed<-renderPlot({
    if(ready()){
      select<-input$exon_table_rows_selected
      p<-plot_exons(tx_df= txdf()$collapsed, select=select)
      p
    } else {
      NULL
    }
  }, height = function() {
    if(ready()){
      50
    } else {
      1
    }
  })
  
  
  ####### junction related ########
  
  draw_junction<-reactive({
    junction_data<-expression_data()$junction
    if(ready() & length(unique(junction_data$junction_id))>0){
      return(T)
    } else {
      return(F)
    }
  })
  
  output$junc_table<-renderDataTable({
    if(draw_junction() & ready()) {
      dat<-unique(expression_data()$junction[,c("junction_id", "start", "end")])
      rownames(dat)<-c(1:nrow(dat))
      datatable(dat)
    } else {
      NULL
    }
  })
  
  output$deselect_junc<-renderUI({
    if(ready() & draw_junction()){
      actionButton(session$ns("deselect_junc_bttn"), label = "Deselect")
    } else {
      NULL
    }
  })
  
  junc_select_proxy<-dataTableProxy("junc_table", session = session)
  observeEvent(input$deselect_junc_bttn, {
    selectRows(junc_select_proxy, NULL)
  })
  
  selected_junctions<-reactive({
    junctions<-unique(expression_data()$junction[,c("junction_id")])
    junctions<-junctions[input$junc_table_rows_selected]
    return(junctions)
  })
  
  output$junc_alert<-renderUI({
    if(!ready()){
      msg<-"Please select genes and tissues above to see expression distribution"
      color<-"maroon" 
      condition<-T
    } else if (ready() && length(selected_junctions())<1 && draw_junction()){
      msg<-"Select exons from the table above to see their expression"
      color<-"light-blue"
      condition<-T
    } else if (ready() && !draw_junction()){
      msg<-"There are no junctions to draw"
      color<-"olive"
      condition<-T
    } else {
      condition<-F
    }
    tagList(
      alert_box(alert = msg, color = color, condition = condition)
    )
  })
  
  output$junc_exp<-renderPlotly({
    if(draw_junction() & ready()){
      plot_data<-expression_data()$junction[(expression_data()$junction$junction_id %in% selected_junctions()), ]
      plot_ly(data = plot_data,
              x=~tissue, y=~value, 
              color=~junction_id, type = "box", source = "junction_boxplot")%>%
        layout(boxmode = "group", xaxis=list(title="Tissue"), 
               yaxis=list(title="Read Count"))
    } else {
      NULL
    }
  })
  
  output$junc_gene_model<-renderPlot({
    if(draw_junction() & ready()) {
      junc_data<-expression_data()$junction
      select<-input$junc_table_rows_selected
      p<-plot_junctions(junc_exp=junc_data, tx_df= txdf()$collapsed, select=select)
      p
    } else {
      NULL
    }
  },height = function() {
    if(ready() && draw_junction()){
      200
    } else {
      1
    }
  })
  
}






