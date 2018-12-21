########## FUNCTIONS ###########

get_exon_df<-function(conn, gene_name=NULL, collapsed=F){
  if(!is.null(gene_name)){
    if(collapsed){
      table<-"annotation.collapsed"
    } else {
      table<-"annotation.full"
    }
    query<-paste0("select * from ", table, " where geneid=?gene")
    sql<-sqlInterpolate(conn, query, gene=gene_name)
    df<-dbGetQuery(conn, sql, gene=gene_name)[,c("geneid", "txname", "exonname", "start", "end", "strand", "exonrank")]
    df$tx_rank<-as.integer(as.factor(df$txname))
    df$class<-"exon"
    actual_start<-min(df$start)
    exons_data<-list(df=df, start=actual_start)
    return(exons_data)
  } else {
    NULL
  }
}

#this needs to be normalized by the gene start
get_introns<-function(df, actual_start, txname){
  tx<-df[which(df$txname==txname), ]
  if(nrow(tx)>2){
    tx<-tx[base::order(tx$start),]
    intron_start<-tx$start[-nrow(tx)]+1
    intron_end<-tx$end[-1]-1
    intron_name<-paste(tx$exonname[-nrow(tx)], tx$exonname[-1], sep="_") 
    len<-length(intron_name)
    introndf<-data.frame(geneid=rep(tx$gene[1], len), exonname=intron_name, 
                         strand=rep(tx$strand[1], len), start=intron_start, 
                         end=intron_end, exonrank=seq(from=1, to=len), 
                         txname=rep(tx$txname[1], len), tx_rank=rep(tx$tx_rank[1], len), 
                         class=rep("intron", len))
    dat<-rbind(tx, introndf)
    dat<-dat[base::order(dat$start), ]
    intron_loc<-which(dat$class=="intron")
  } else {
    dat<-tx
  }
  return(dat)
}

get_introns<-Vectorize(get_introns, vectorize.args = "txname", SIMPLIFY = F)

parse_exon_df<-function(conn, gene_name, collapsed){
  if(!is.null(gene_name)){
    df<-get_exon_df(conn = conn, gene_name = gene_name, collapsed = collapsed)
    txs<-unique(df$df$txname)
    tx_df<-do.call("rbind", get_introns(df$df, txname = txs, actual_start = df$start))
    return(tx_df)
  }else {
    NULL
  }
}

plot_transcripts<-function(tx_df, hover=NULL){
  if(is.null(hover)){
    p<-ggplot()+
      geom_rect(data=tx_df[tx_df$class=="exon",],
                aes(xmin = start, xmax = end, ymax = tx_rank + 0.25, 
                    ymin = -tx_rank - 0.25))+
      geom_rect(data=tx_df[tx_df$class=="intron",],
                aes(xmin = start, xmax = end, ymax = -tx_rank + 0.01, 
                    ymin = tx_rank - 0.01))+
      theme_void()+guides(fill=F, color=F)
  } else {
    p<-ggplot()+
      geom_rect(data=tx_df[tx_df$class=="exon",],
                aes(xmin = start, xmax = end, ymax = tx_rank + 0.25, 
                    ymin = tx_rank - 0.25), alpha=0.3)+
      geom_rect(data=tx_df[tx_df$class=="intron",],
                aes(xmin = start, xmax = end, ymax = tx_rank + 0.01, 
                    ymin = tx_rank - 0.01), alpha=0.3)+
      # new ggproto on top
      geom_rect(data=tx_df[tx_df$class=="exon" & tx_df$tx_rank==hover,],
                aes(xmin = start, xmax = end, ymax = tx_rank + 0.25, 
                    ymin = tx_rank - 0.25), fill="firebrick")+
      geom_rect(data=tx_df[tx_df$class=="intron"& tx_df$tx_rank==hover,],
                aes(xmin = start, xmax = end, ymax = tx_rank + 0.01, 
                    ymin = tx_rank - 0.01), color="firebrick")+
      theme_void()+guides(fill=F, color=F)
  }
  return(p)
}

# need to add a function that links the plotly to the ggplot for brush events
plot_exons<-function(tx_df, select){
  if(is.null(select)){
    p<-ggplot()+
      geom_rect(data=tx_df[tx_df$class=="exon",],
                aes(xmin = start, xmax = end, ymax = tx_rank + 0.25, 
                    ymin = tx_rank - 0.25))+
      geom_rect(data=tx_df[tx_df$class=="intron",],
                aes(xmin = start, xmax = end, ymax = tx_rank + 0.01, 
                    ymin = tx_rank - 0.01))+
      theme_void()+guides(fill=F, color=F)
  } else {
    exondf<-tx_df[tx_df$class=="exon",]
    p<-ggplot()+
      geom_rect(data=tx_df[tx_df$class=="exon",],
                aes(xmin = start, xmax = end, ymax = tx_rank + 0.25, 
                    ymin = tx_rank - 0.25), alpha=0.3)+
      geom_rect(data=tx_df[tx_df$class=="intron",],
                aes(xmin = start, xmax = end, ymax = tx_rank + 0.01, 
                    ymin = tx_rank - 0.01), alpha=0.3)+
      # new ggproto on top
      geom_rect(data=exondf[select,],
                aes(xmin = start, xmax = end, ymax = tx_rank + 0.25, 
                    ymin = tx_rank - 0.25), fill="firebrick")+
      theme_void()+guides(fill=F, color=F)
  }
  return(p)
}

plot_junctions<-function(junc_exp, tx_df, select){
  uniq_junc<-unique(junc_exp[, 1:3])
  uniq_junc$start<-as.integer(uniq_junc$start)
  uniq_junc$end<-as.integer(uniq_junc$end)
  if(length(select)==0){
    p<-ggplot()+
      geom_rect(data=tx_df[tx_df$class=="exon",],
                aes(xmin = start, xmax = end, ymax = 0.1, 
                    ymin = 0), alpha=0.3)+
      geom_curve(data=uniq_junc, aes(x=start, xend=end, y=0.1, yend=0.1), 
                 curvature=-0.3)+ylim(0, 1)+
      theme_void()+guides(fill=F, color=F)
  } else {
    p<-ggplot()+
      geom_rect(data=tx_df[tx_df$class=="exon",],
                aes(xmin = start, xmax = end, ymax = 0.1, 
                    ymin = 0), alpha=0.3)+
      geom_curve(data=uniq_junc, aes(x=start, xend=end, y=0.1, yend=0.1), 
                 curvature=-0.3, alpha=0.3)+
      geom_curve(data=uniq_junc[select,], 
                 aes(x=start, xend=end, y=0.1, yend=0.1), 
                 curvature=-0.3, color="firebrick", size=3)+ylim(0, 1)+
      theme_void()+guides(fill=F, color=F)
  }
  return(p)
}


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
    dat<-txdf_collapsed()[txdf_collapsed()$class=="exon",]
    dat<-dat[order(as.integer(dat$exonrank)),]
    rownames(dat)<-c(1:nrow(dat))
    datatable(dat[,c("exonname", "exonrank")])
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





