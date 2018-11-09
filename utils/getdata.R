# so far this works for gene and isoform
get_expression<-function(datadir, gene_id, tissue, samples, table, extra_columns=NULL){
  tiss_samples<-samples[[tissue]]$SAMPID
  columns<-c("gene_id", tiss_samples)
  if(!is.null(extra_columns)){
    columns<-c(extra_columns, columns)
  }
  columns<-paste0("[", columns, "]", collapse = ",")
  expressions<-list()
  #for(tissue in tissues){
  name<-stripWhitespace(tissue)
  name<-gsub(" ", "_", name)
  db_loc<-paste0(datadir, "databases/", name, "_gtex.db")
  db<-dbConnect(drv = SQLite(), dbname=db_loc)
  if(length(gene_id)>1){
    gene_id<-as.character(gene_id)
    genes_combined<-paste("'", gene_id, "'", sep="", collapse = ",")
    query<-paste("select ", columns, " from ", table,  " where [gene_id] in (", 
                 genes_combined, ")", sep="")
  } else {
    query<-paste0("select ", columns, " from ", table, " where gene_id='", gene_id, "'")
  }
  results<-dbGetQuery(db, query)
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









