median_heatmap_ui<-function(id){
  ns<-NS(id)
  plotlyOutput(ns("tpm_heat"), height = "600px")
}

median_heatmap<-function(session, input, output, genes_df){
  if(!(is.null(gene_df()))){
    forheat<-genes_df()$genes
    tissues<-colnames(forheat)[-c(1:3)]
    genes<-forheat$genesymbol
    mat<-as.matrix(forheat[,-c(1:3)])
    ## rescale colors
    vals <- unique(scales::rescale(c(mat)))
    o <- order(vals, decreasing = FALSE)
    cols <- scales::col_numeric("Spectral", domain = NULL)(vals)
    colz <- setNames(data.frame(vals[o], cols[o]), NULL)
    ####
    plot_ly(y=tissues, x=genes, z=t(mat), type="heatmap", 
            source="tpm_heatplot", colorscale=colz) %>% 
      layout(xaxis=list(title="", dtick=1), 
             yaxis=list(title="", dtick=1))
  } else {
    NULL
  }
}