gene_expression_ui<-function(id){
  ns<-NS(id)
  uiOutput(ns("gene_expression_plot"))
}

#this needs to be paginated
gene_expression<-function(input, output, session, login, conn, genes, tissues_samples){
  
  
}