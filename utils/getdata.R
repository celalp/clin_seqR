# so far this works for gene and isoform
get_expression<-function(conn, gene_id, tissue, samples, table, extra_columns=NULL){
  tiss_samples<-tolower(samples[[tissue]]$sampid)
  tiss_samples<-paste('"', tiss_samples, '"', sep = "")
  columns<-c("gene_id", tiss_samples)
  if(!is.null(extra_columns)){
    columns<-c(extra_columns, columns)
  }
  columns<-paste0(columns, collapse = ",")
  expressions<-list()
  #for(tissue in tissues){
  tissue<-gsub(" ", "_", tissue)
  table<-paste(tissue, table, sep=".")
  if(length(gene_id)>1){
    gene_id<-as.character(gene_id)
    genes_combined<-paste("'", gene_id, "'", sep="", collapse = ",")
    query<-paste("select ", columns, " from ", table,  " where gene_id in (", 
                 genes_combined, ")", sep="")
    print(query)
  } else {
    query<-paste0("select ", columns, " from ", table, " where gene_id='", gene_id, "'")
  }
  results<-dbGetQuery(conn, query)
  if(nrow(results)==0){ #this is for junctions
    results<-NULL
    expressions[[tissue]]<-results
  } else {
    results$tissue<-tissue
    expressions[[tissue]]<-melt(results)
  }
  return(expressions)
}

get_expression<-Vectorize(get_expression, vectorize.args = c("tissue"), SIMPLIFY = T)









