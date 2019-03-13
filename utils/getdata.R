
# so far this works for gene and isoform
get_expression<-function(conn, gene_id, tissues_samples, table, extra_columns=NULL){
  tissues<-names(tissues_samples)
  expressions<-list()
  for(tiss in tissues){
    tiss_samples<-tolower(tissues_samples[[tiss]])
    tiss_samples<-paste('"', tiss_samples, '"', sep = "")
    columns<-c("gene_id", tiss_samples)
    if(!is.null(extra_columns)){
      columns<-c(extra_columns, columns)
    }
    columns<-paste0(columns, collapse = ",")
    tablename<-paste(tiss, table, sep=".")
    if(length(gene_id)>1){
      gene_id<-as.character(gene_id)
      genes_combined<-paste("'", gene_id, "'", sep="", collapse = ",")
      query<-paste("select ", columns, " from ", tablename,  " where gene_id in (", 
                   genes_combined, ")", sep="")
    } else {
      query<-paste0("select ", columns, " from ", tablename, " where gene_id='", gene_id, "'")
    }
    results<-dbGetQuery(conn, query)
    if(nrow(results)==0){ #this is for junctions
      results<-NULL
      expressions[[tiss]]<-results
    } else {
      results$tissue<-tiss
      expressions[[tiss]]<-melt(results)
    }
  }
  expressions<-do.call("rbind", expressions)
  return(expressions)
}

#get_expression<-Vectorize(get_expression, vectorize.args = c("tissue"), SIMPLIFY = T)

# filter functions

get_samples_subjects<-function(conn, tissues){
  if(!is.null(tissues)){
    tissue_select<-paste0("'", tissues, sep = "'")
    tissue_select<-paste(tissue_select, collapse = ",")
    query<-paste0("select * from samples.samples where smtsd in (", tissue_select, ")")
    samples_subjects<-dbGetQuery(conn, query)
    samples_subjects$sex<-ifelse(samples_subjects$sex==1, "Male", "Female")
    samples_subjects$dthhrdy<-as.character(samples_subjects$dthhrdy)
    samples_subjects$dthhrdy<-gsub("0", "Ventilator Case", samples_subjects$dthhrdy)
    samples_subjects$dthhrdy<-gsub("1", "Fast and Violent", samples_subjects$dthhrdy)
    samples_subjects$dthhrdy<-gsub("2", "Fast", samples_subjects$dthhrdy)
    samples_subjects$dthhrdy<-gsub("3", "Intermediate", samples_subjects$dthhrdy)
    samples_subjects$dthhrdy<-gsub("4", "Slow", samples_subjects$dthhrdy)
    samples_subjects$dthhrdy[which(is.na(samples_subjects$dthhrdy))]<-"Unavailable"
    return(samples_subjects)
  } else {
    NULL
  }
}

get_common<-function(samples_subjects){
  each_tissue<-split(samples_subjects[,c("smtsd", "subjid")], as.factor(samples_subjects$smtsd))
  common<-plyr::join_all(each_tissue, by="subjid", type="inner")$subjid
}

# genevis functions

get_exon_df<-function(conn, gene_name=NULL, collapsed=F){
  if(!is.null(gene_name)){
    if(collapsed){
      table<-"annotation.gtex_collapsed"
    } else {
      table<-"annotation.gtex_full"
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

# trying to be dry need tests for this stuff
get_tissues_samples<-function(table){
  samples<-table$sampid
  tissues<-table$smtsd
  tissues<-gsub("\\(", " ", tissues)
  tissues<-gsub("\\)", " ", tissues)
  tissues<-gsub("-", " ", tissues)
  tissues<-stripWhitespace(tissues)
  tissues<-gsub(" $", "", tissues)
  tissues<-gsub(" ", "_", tissues)
  tissues<-tolower(tissues)
  sample_list<-split(samples, as.factor(tissues))
  return(sample_list)
}


