median_heatmap_ui<-function(id){
  ns<-NS(id)
  uiOutput(ns("tpm_heat_ui"), height = "600px")
}

median_heatmap<-function(input, output, session, conn, login, genesymbol){
  
  output$tpm_heat_ui<-renderUI({
    if(login$login){
      if(!is.null(genesymbol()$genesymbol) | length(genesymbol()$genesymbol >0 )){
        genes<-paste("'", genesymbol()$genesymbol, "'", sep="", collapse = ",")
        query<-paste("select * from samples.median_expression where genesymbol in (", genes, ")", sep = "")
        query<-sqlInterpolate(conn, query)
        forheat<-dbGetQuery(conn, query)
        tissues<-colnames(forheat)[-c(1:3)]
        genes<-as.character(forheat$genesymbol)
        mat<-as.matrix(forheat[,-c(1:3)])
        ## rescale colors
        vals <- unique(scales::rescale(c(mat)))
        o <- order(vals, decreasing = FALSE)
        cols <- scales::col_numeric("Spectral", domain = NULL)(vals)
        colz <- setNames(data.frame(vals[o], cols[o]), NULL)
        ####
        output$median_heatmap<-renderPlotly({
          plot_ly(y=tissues, x=genes, z=t(mat), type="heatmap", 
                  source="tpm_heatplot", colorscale=colz) %>% 
            layout(xaxis=list(title="", dtick=1), 
                   yaxis=list(title="", dtick=1))
        })
        tagList(
          plotlyOutput(session$ns("median_heatmap"), height = "600px")
        )
      } else {
        #TODO notification
        NULL
      }
    } else {
      NULL
    }
  })
}