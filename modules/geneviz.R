
genevis_ui<-function(id){
  ns<-NS(id)
  
  tagList(
    uiOutput(ns("gene_picker_ui")),
    tags$br(),
    tabsetPanel( # need to mess with the colors of the tab labels
      tabPanel("Isoform Expression", value = "gtex_isoform",
               column(width=2, 
                      uiOutput(ns("selection_ui"))),
               column(width=10,
                      tagList(
                        plotlyOutput(ns("isoform_exp")),
                        plotOutput(ns("isoform_model"))
                      )
               )
      ),
      tabPanel("Exon Expression", value = "gtex_exon", 
               fluidRow(
                 column(width=4,
                        tags$br(),
                        dataTableOutput(ns("exon_table")),
                        actionButton(ns("deselect_exon"), label = "Deselect")
                 ),
                 column(width=7, offset = 1, 
                        tags$br(),
                        plotlyOutput(ns("exon_exp"))
                 )),
               tags$hr(),
               plotOutput(ns("gene_model_collapsed"))
      ),
      tabPanel("Junction Expression", value="gtex_junc", 
               fluidRow(
                 column(width=4, 
                        tags$br(), 
                        dataTableOutput(ns("junc_table")),
                        actionButton(ns("deselect_junc"), label = "Deselect")
                 ),
                 column(width=7, offset = 1, 
                        tags$br(),
                        plotlyOutput(ns("junc_exp"))
                 )
               ),
               tags$hr(),
               plotOutput(ns("junc_gene_model"))
      )
    )
  )
}

genevis<-function(input, output, session, tissues_samples, genes, login, conn){ # this needs data as well
  
  output$gene_picker_ui<-renderUI({
    if(login$login){
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
    if(login$login){
      tagList(
        tags$br(),
        radioGroupButtons(session$ns("isoform_tpm_reads"), label = "TPM or read count?", 
                          choices = c("TPM", "Read Count"), direction = "vertical", 
                          selected = "TPM"))
    } else {
      NULL
    }
  })
  
  txdf<-reactive({
    gene_name<-genes()$geneid[genes()$genesymbol==input$gene_picker]
    txdf_full<-parse_exon_df(conn = conn, gene_name = gene_name, collapsed=F)
    txdf_collapsed<-parse_exon_df(conn = conn, gene_name = gene_name, collapsed=T)
    txdf<-list(gene_id=gene_name, full=txdf_full, collapsed=txdf_collapsed)
    return(txdf)
  })
  
  expression_data<-reactive({
    if(!is.null(txdf()$gene_id)){
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
      junction_expression<<-get_expression(conn=conn, gene_id=txdf()$gene_id, 
                                          tissues_samples=tissues_samples(), 
                                          table="junction_reads", extra_columns=c('"junction_id"', '"start"','"end"'))
      expression_data<-list(isoform=isoform_expression, exon=exon_expression, junction=junction_expression)
    } else {
      expression_data<-list(isoform=NULL, exon=NULL, junction=NULL)
    }
    return(expression_data)
  })
  
  
  output$isoform_exp<-renderPlotly({
    if(login$login & length(txdf()$gene_id)>0){
        plot_ly(data = expression_data()$isoform, x=~tissue, y=~value, 
                color=~transcript_id, type = "box", source = "isoform_boxplot")%>%
          layout(boxmode = "group", xaxis=list(title="Tissue"), 
                 yaxis=list(title=input$isoform_tpm_reads))
      } else {
        NULL
      }
  })
  
  output$isoform_model<-renderPlot({
    if(login$login & length(txdf()$gene_id)>0){
      hover<-event_data(event = "plotly_hover", source="isoform_boxplot")$curveNumber+1
      p<-plot_transcripts(tx_df= txdf()$full, hover=hover)
      p
    } else {
      NULL
    }
  },height = function() {
    if(nrow(txdf()$full)>0){
      max(txdf()$full$tx_rank*40)
    } else {
      50
    }
  })
  
  output$exon_table<-renderDataTable({
    if(login$login & length(txdf()$gene_id)>0){
      dat<-txdf()$collapsed
      dat<-dat[dat$class=="exon",]
      dat<-dat[order(as.integer(dat$exonrank)),]
      rownames(dat)<-c(1:nrow(dat))
      datatable(dat[,c("exonname", "exonrank")])
    } else {
      NULL
    }
  })
  
  #TODO make this a ui
  exon_select_proxy<-dataTableProxy("exon_table", session = session)
  observeEvent(input$deselect_exon, {
    selectRows(exon_select_proxy, NULL)
  })
  
  output$exon_exp<-renderPlotly({
    nointron<-txdf()$collapsed[txdf()$collapsed$class=="exon",]
    exons<-nointron$exonname[input$exon_table_rows_selected]
    if(login$login & length(txdf()$gene_id)>0){
      if (length(exons)>0){
        plot_data<-expression_data()$exon[(expression_data()$exon$exon_id %in% exons), ]
        plot_ly(data = plot_data,
                x=~tissue, y=~value, 
                color=~exon_id, type = "box", source = "exon_boxplot")%>%
          layout(boxmode = "group", xaxis=list(title="Tissue"), 
                 yaxis=list(title="Read Count"))
      } else {
        NULL
      }
    }else {
      NULL
    }
  })
  
  #need to change the highlighting mode
  output$gene_model_collapsed<-renderPlot({
    if(login$login & length(txdf()$gene_id)>0){
      select<-input$exon_table_rows_selected
      p<-plot_exons(tx_df= txdf()$collapsed, select=select)
      p
    } else {
      NULL
    }
  }, height = 50)
  
  
  draw_junction<-reactive({
    if(!is.null(txdf()$gene_id) && dim(do.call("rbind",expression_data()$junction))[2]>0){
      return(T)
    } else {
      return(F)
    }
  })
  
  
  output$junc_table<-renderDataTable({
    dat<-unique(expression_data()$junction[,c("junction_id", "start", "end")])
    rownames(dat)<-c(1:nrow(dat))
    datatable(dat)
  })
  
  junc_select_proxy<-dataTableProxy("junc_table", session = session)
  observeEvent(input$deselect_junc, {
    selectRows(junc_select_proxy, NULL)
  })
  
  
  output$junc_exp<-renderPlotly({
    junctions<-unique(expression_data()$junction[,c("junction_id", "start", "end")])
    selected<-junctions$junction_id[input$junc_table_rows_selected]
    if(draw_junction()){
      if (length(selected)>0){
        plot_data<-expression_data()$junction[(expression_data()$junction$junction_id %in% selected), ]
        plot_ly(data = plot_data,
                x=~tissue, y=~value, 
                color=~junction_id, type = "box", source = "junction_boxplot")%>%
          layout(boxmode = "group", xaxis=list(title="Tissue"), 
                 yaxis=list(title="Read Count"))
      } else {
        NULL
      }
    }else {
      NULL
    }
  })
  
  
  #need to change the highlighting mode
  output$junc_gene_model<-renderPlot({
    if(draw_junction()) {
      junc_data<-expression_data()$junction
      select<-input$junc_table_rows_selected
      p<-plot_junctions(junc_exp=junc_data, tx_df= txdf()$collapsed, select=select)
      p
    } else {
      NULL
    }
  }, height = 200)
  
}






