########## FUNCTIONS ###########




########## MODULE ###########

genevis_ui<-function(id){
  ns<-NS(id)
  tabsetPanel( # need to mess with the colors of the tab labels
    tabPanel("Isoform Expression", value = "gtex_isoform",
             column(width=2, 
                    radioGroupButtons(ns("isoform_tpm_reads"), label = "TPM or read count?", 
                                      choices = c("TPM", "Read Count"), direction = "vertical", 
                                      selected = "TPM")),
             column(width=10,
                    #withSpinner(
                    tagList(
                      plotlyOutput(ns("isoform_exp")),
                      plotOutput(ns("gene_model_exp"))
                      #  color = '#202020'
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
                      bsAlert("exon_alert"),
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
                      bsAlert("junc_alert"),
                      plotlyOutput(ns("junc_exp"))
               )
             ),
             tags$hr(),
             plotOutput(ns("junc_gene_model"))
    )
  )
}

genevis<-function(input, output, session, datadir, tissues, samples, gene_id,
                  conn){ # this needs data as well
  
  txdf_full<-reactive({
    txdf<-parse_exon_df(conn = conn, gene_name = gene_id()[1], collapsed=F)
    return(txdf)
  })
  
  txdf_collapsed<-reactive({
    txdf<-parse_exon_df(conn = conn, gene_name = gene_id()[1], collapsed = T)
    return(txdf)
  })
  
  
  expression_data<-reactive({
    if(!is.null(gene_id())){
      if(input$isoform_tpm_reads=="TPM"){
        iso_tbl<-"transcript_tpm"
      } else {
        iso_tbl<-"transcript_reads"
      }
      isoform_expression<-get_expression(conn=conn, gene_id=gene_id()[1], 
                                         samples=samples(), tissue=tissues(), 
                                         table=iso_tbl, extra_columns='"transcript_id"')
      isoform_expression<-do.call("rbind", isoform_expression)
      exon_expression<-get_expression(conn=conn, gene_id=gene_id()[1], 
                                      samples=samples(), tissue=tissues(), 
                                      table="exon_reads", extra_columns='"exon_id"')
      exon_expression<-do.call("rbind", exon_expression)
      junction_expression<-get_expression(conn=conn, gene_id=gene_id()[1], 
                                          samples=samples(), tissue=tissues(), 
                                          table="junctions", extra_columns=c('"junction_id"', '"start"','"end"'))
      junction_expression<-do.call("rbind", junction_expression)
      expression_data<-list(isoform=isoform_expression, exon=exon_expression, junction=junction_expression)
    } else {
      expression_data<-list(isoform=NULL, exon=NULL, junction=NULL)
    }
    return(expression_data)
  })
  
  output$isoform_exp<-renderPlotly({
    if(!is.null(gene_id())){
      plot_ly(data = expression_data()$isoform, x=~tissue, y=~value, 
              color=~transcript_id, type = "box", source = "isoform_boxplot")%>%
        layout(boxmode = "group", xaxis=list(title="Tissue"), 
               yaxis=list(title=input$isoform_tpm_reads))
    } else {
      NULL
    }
  })
  
  output$gene_model_exp<-renderPlot({
    if(!is.null(gene_id())){
      hover<-event_data(event = "plotly_hover", source="isoform_boxplot")$curveNumber+1
      p<-plot_transcripts(tx_df= txdf_full(), hover=hover)
      p
    } else {
      NULL
    }
  },height = function() {
    if(!is.null(txdf_full())){
      max(txdf_full()$tx_rank*40)
    } else {
      50
    }
  }
  )
  
  output$exon_table<-renderDataTable({
    if(is.null(gene_id())){
      NULL
    } else {
    dat<-txdf_collapsed()[txdf_collapsed()$class=="exon",]
    dat<-dat[order(as.integer(dat$exonrank)),]
    rownames(dat)<-c(1:nrow(dat))
    datatable(dat[,c("exonname", "exonrank")])
    }
  })
  
  exon_select_proxy<-dataTableProxy("exon_table", session = session)
  observeEvent(input$deselect_exon, {
    selectRows(exon_select_proxy, NULL)
  })
  
  output$exon_exp<-renderPlotly({
    nointron<-txdf_collapsed()[txdf_collapsed()$class=="exon",]
    exons<-nointron$exonname[input$exon_table_rows_selected]
    if(!is.null(gene_id())){
      if (length(exons)>0){
        closeAlert(session, "exon_alert_control")
        plot_data<-expression_data()$exon[(expression_data()$exon$exon_id %in% exons), ]
        plot_ly(data = plot_data,
                x=~tissue, y=~value, 
                color=~exon_id, type = "box", source = "exon_boxplot")%>%
          layout(boxmode = "group", xaxis=list(title="Tissue"), 
                 yaxis=list(title="Read Count"))
      } else {
        createAlert(session, "exon_alert", "exon_alert_control", title = "",
                    content = "Select exons on the left to see their expression", append = FALSE)
      }
    }else {
      NULL
    }
  })
  
  #need to change the highlighting mode
  output$gene_model_collapsed<-renderPlot({
    if(!is.null(gene_id())){
      select<-input$exon_table_rows_selected
      p<-plot_exons(tx_df= txdf_collapsed(), select=select)
      p
    } else {
      NULL
    }
  }, height = 50)
  
  draw_junction<-reactive({
    if(!is.null(gene_id()) && dim(do.call("rbind",expression_data()$junction))[2]>0){
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
        closeAlert(session, "junction_alert_control")
        plot_data<-expression_data()$junction[(expression_data()$junction$junction_id %in% selected), ]
        plot_ly(data = plot_data,
                x=~tissue, y=~value, 
                color=~junction_id, type = "box", source = "junction_boxplot")%>%
          layout(boxmode = "group", xaxis=list(title="Tissue"), 
                 yaxis=list(title="Read Count"))
      } else {
        createAlert(session, "junc_alert", "junction_alert_control", title = "",
                    content = "Select junction on the left to see their expression", append = FALSE)
      }
    }else {
      toastr_warning("There is no junction data associated with this gene")
    }
  })

  
  #need to change the highlighting mode
  output$junc_gene_model<-renderPlot({
    if(draw_junction()) {
      junc_data<-expression_data()$junction
      select<-input$junc_table_rows_selected
      p<-plot_junctions(junc_exp=junc_data, tx_df= txdf_collapsed(), select=select)
      p
    } else {
      NULL
    }
  }, height = 200)
}





