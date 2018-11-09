########## FUNCTIONS ###########

get_exon_df<-function(annot, gene_name){
  df<-AnnotationDbi::select(annot, keys=gene_name, 
                            columns=c("GENEID", "TXNAME", "EXONNAME", "EXONRANK", 
                                      "EXONSTART", "EXONEND", "EXONSTRAND"),
                            keytype = "GENEID")
  df$rank<-as.integer(as.factor(df$TXNAME))
  colnames(df)<-c("gene", "name", "strand", "start", "end", "rank", "tx", "tx_rank")
  df$class<-"exon"
  actual_start<-min(df$start)
  exons_data<-list(df=df, start=actual_start)
  return(exons_data)
}

#this needs to be normalized by the gene start
get_introns<-function(df, actual_start, txname){
  tx<-df[which(df$tx==txname), ]
  if(nrow(tx)>2){
    tx<-tx[base::order(tx$start),]
    intron_start<-tx$start[-nrow(tx)]+1
    intron_end<-tx$end[-1]-1
    intron_name<-paste(tx$name[-nrow(tx)], tx$name[-1], sep="_") 
    len<-length(intron_name)
    introndf<-data.frame(gene=rep(tx$gene[1], len), name=intron_name, 
                         strand=rep(tx$strand[1], len), start=intron_start, 
                         end=intron_end, rank=seq(from=1, to=len), 
                         tx=rep(tx$tx[1], len), tx_rank=rep(tx$tx_rank[1], len), 
                         class=rep("intron", len))
    dat<-rbind(tx, introndf)
    dat<-dat[base::order(dat$start), ]
    intron_loc<-which(dat$class=="intron")
  } else {
    dat<-tx
  }
  lengths<-dat$end-dat$start
  relative_start<-min(dat$start)-actual_start
  dat$new_start<-c(relative_start, (relative_start + cumsum(lengths[-length(lengths)])))
  dat$new_end<-cumsum(lengths)+relative_start
  return(dat)
}

get_introns<-Vectorize(get_introns, vectorize.args = "txname", SIMPLIFY = F)

parse_exon_df<-function(annot, gene_name){
  df<-get_exon_df(annot = annot, gene_name = gene_name)
  txs<-unique(df$df$tx)
  tx_df<-do.call("rbind", get_introns(df$df, txname = txs, actual_start = df$start))
  return(tx_df)
}

plot_transcripts<-function(tx_df, hover=NULL){
  if(is.null(hover)){
    p<-ggplot()+
      geom_rect(data=tx_df[tx_df$class=="exon",],
                aes(xmin = new_start, xmax = new_end, ymax = tx_rank + 0.25, 
                    ymin = -tx_rank - 0.25, text=tx))+
      geom_rect(data=tx_df[tx_df$class=="intron",],
                aes(xmin = new_start, xmax = new_end, ymax = -tx_rank + 0.01, 
                    ymin = tx_rank - 0.01, text=tx))+
      theme_void()+guides(fill=F, color=F)
  } else {
    p<-ggplot()+
      geom_rect(data=tx_df[tx_df$class=="exon",],
                aes(xmin = new_start, xmax = new_end, ymax = tx_rank + 0.25, 
                    ymin = tx_rank - 0.25), alpha=0.3)+
      geom_rect(data=tx_df[tx_df$class=="intron",],
                aes(xmin = new_start, xmax = new_end, ymax = tx_rank + 0.01, 
                    ymin = tx_rank - 0.01), alpha=0.3)+
      # new ggproto on top
      geom_rect(data=tx_df[tx_df$class=="exon" & tx_df$tx_rank==hover,],
                aes(xmin = new_start, xmax = new_end, ymax = tx_rank + 0.25, 
                    ymin = tx_rank - 0.25), fill="firebrick")+
      geom_rect(data=tx_df[tx_df$class=="intron"& tx_df$tx_rank==hover,],
                aes(xmin = new_start, xmax = new_end, ymax = tx_rank + 0.01, 
                    ymin = tx_rank - 0.01), color="firebrick")+
      theme_void()+guides(fill=F, color=F)
  }
  return(p)
}

# need to add a function that links the plotly to the ggplot for brush events
plot_exons<-function(tx_df, hover, select=NULL){
  if(is.null(select)){
    #tx_df<-txdf
  } else {
    max_curve<-max(select$curveNumber)
    min_curve<-min(select$curveNumber)
  }
  if(is.null(hover)){
    p<-ggplot()+
      geom_rect(data=tx_df[tx_df$class=="exon",],
                aes(xmin = new_start, xmax = new_end, ymax = tx_rank + 0.25, 
                    ymin = tx_rank - 0.25, text=tx))+
      geom_rect(data=tx_df[tx_df$class=="intron",],
                aes(xmin = new_start, xmax = new_end, ymax = tx_rank + 0.01, 
                    ymin = tx_rank - 0.01, text=tx))+
      theme_void()+guides(fill=F, color=F)
  } else {
    p<-ggplot()+
      geom_rect(data=tx_df[tx_df$class=="exon",],
                aes(xmin = new_start, xmax = new_end, ymax = tx_rank + 0.25, 
                    ymin = tx_rank - 0.25), alpha=0.3)+
      geom_rect(data=tx_df[tx_df$class=="intron",],
                aes(xmin = new_start, xmax = new_end, ymax = tx_rank + 0.01, 
                    ymin = tx_rank - 0.01), alpha=0.3)+
      # new ggproto on top
      geom_rect(data=tx_df[tx_df$class=="exon" & tx_df$rank==hover,],
                aes(xmin = new_start, xmax = new_end, ymax = tx_rank + 0.25, 
                    ymin = tx_rank - 0.25), fill="firebrick")+
      theme_void()+guides(fill=F, color=F)
  }
  return(p)
}

plot_junctions<-function(junc_exp, tx_df, hover){
  uniq_junc<-unique(junc_exp[, 1:3])
  uniq_junc$start<-as.integer(uniq_junc$start)
  uniq_junc$end<-as.integer(uniq_junc$end)
  if(length(hover)==0){
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
      geom_curve(data=uniq_junc[hover,], 
                 aes(x=start, xend=end, y=0.1, yend=0.1), 
                 curvature=-0.3, color="firebrick", size=3)+ylim(0, 1)+
      theme_void()+guides(fill=F, color=F)
  }
  return(p)
}

#this is for novel isoforms and stuff
get_altered_coord<-function(df, coord, shortened=T, shortened_by=10){
  distances<-coord-df$start
  prec<-distances[distances>0] # find the preceeding ones
  closest<-min(prec)
  loc<-df[which(distances==closest),] # this is the index 
  ratio<-(loc$end-loc$start)/(coord-loc$start)
  new_coord<-loc$new_start+(loc$new_end-loc$new_start)*ratio
  return(new_coord)
}

get_altered_coord<-Vectorize(get_altered_coord, vectorize.args = "coord")


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
             #withSpinner(
             column(width=10, offset = 1,
                    plotlyOutput(ns("exon_exp")),
                    plotOutput(ns("gene_model_collapsed"))
             )
    ),
    tabPanel("Junction Expression", value="gtex_junc", 
             column(width = 10, offset = 1, 
                    plotlyOutput(ns("junc_exp")),
                    plotOutput(ns("junc_gene_model"))
             )
    )
  )
}

genevis<-function(input, output, session, datadir, tissues, samples,
                  annot, collapsed_annot, gene_id){ # this needs data as well
  
  txdf_transcript<-reactive({
    txdf<-parse_exon_df(annot = annot, gene_name = gene_id()[1,1])
    return(txdf)
  })
  
  txdf_exon<-reactive({
    txdf<-parse_exon_df(annot = collapsed_annot, gene_name = gene_id()[1,1])
    return(txdf)
  })
  
  txdf_junction<-reactive({
    txdf<-parse_exon_df(annot = collapsed_annot, gene_name = gene_id()[1,1])
    return(txdf)
  })
  
  expression_data<-reactive({
    if(!is.null(gene_id())){
      if(input$isoform_tpm_reads=="TPM"){
        iso_tbl<-"transcript_tpm"
      } else {
        iso_tbl<-"transcript_reads"
      }
      isoform_expression<-get_expression(datadir=datadir, gene_id=gene_id()[1], 
                                         samples=samples(), tissue=tissues(), 
                                         table=iso_tbl, extra_columns="transcript_id")
      exon_expression<-get_expression(datadir=datadir, gene_id=gene_id()[1], 
                                      samples=samples(), tissue=tissues(), 
                                      table="exon_reads", extra_columns="exon_id")
      junction_expression<-get_expression(datadir=datadir, gene_id=gene_id()[1], 
                                          samples=samples(), tissue=tissues(), 
                                          table="junctions", extra_columns=c("junction_id", "start","end"))
      expression_data<-list(isoform=isoform_expression, exon=exon_expression, junction=junction_expression)
    } else {
      expression_data<-list(isoform=NULL, exon=NULL, junction=NULL)
    }
    return(expression_data)
  })
  
  output$isoform_exp<-renderPlotly({
    if(!is.null(gene_id())){
      plot_ly(data = do.call("rbind", expression_data()$isoform), x=~tissue, y=~value, 
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
      p<-plot_transcripts(tx_df= txdf_transcript(), hover=hover)
      p
    } else {
      NULL
    }
  })
  
  output$exon_exp<-renderPlotly({
    if(!is.null(gene_id())){
      plot_ly(data = do.call("rbind", expression_data()$exon), x=~tissue, y=~value, 
              color=~exon_id, type = "box", source = "exon_boxplot")%>%
        layout(boxmode = "group", xaxis=list(title="Tissue"), 
               yaxis=list(title="Read Count"))
    } else {
      NULL
    }
  })
  
  output$gene_model_collapsed<-renderPlot({
    if(!is.null(gene_id())){
      hover<-event_data(event = "plotly_hover", source="exon_boxplot")$curveNumber+1
      p<-plot_exons(tx_df= txdf_exon(), hover=hover)
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
  
  output$junc_exp<-renderPlotly({
    if(draw_junction()){
      plot_ly(data = do.call("rbind", expression_data()$junction), x=~tissue, y=~value, 
              color=~junction_id, type = "box", source = "junction_boxplot")%>%
        layout(boxmode = "group", xaxis=list(title="Tissue"), 
               yaxis=list(title="Read Count"))
    } else {
      toastr_warning("There is no junction data associated with this gene")
    }
  })
  
  output$junc_gene_model<-renderPlot({
    if(draw_junction()) {
      junc_data<-do.call("rbind", expression_data()$junction)
      hover<-unique(event_data(event = "plotly_hover", source="junction_boxplot")$curveNumber+1)
      p<-plot_junctions(junc_exp=junc_data, tx_df= txdf_junction(), hover=hover)
      p
    } else {
      NULL
    }
  }, height = 200)
}





