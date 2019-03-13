# TODO there is a lot of if(is.null) here I feel like this can be restructrured to be more DRY
# so far it is working. 


gene_expression_ui<-function(id){
  ns<-NS(id)
  tagList(
    br(),
    column(width=2, 
           uiOutput(ns("read_tpm_buttons"))
    ), 
    column(width=10, 
           uiOutput(ns("gene_expression_plot"))
    )
  )
}

#this needs to be paginated
gene_expression<-function(input, output, session, login, conn, genes, tissues_samples){
  
  output$read_tpm_buttons<-renderUI({
    if(login$login & length(genes()$geneid>0)){
      radioGroupButtons(session$ns("gene_tpm_reads"), label = "TPM or read count?", 
                        choices = c("TPM", "Read Count"), selected = "TPM",
                        direction = "vertical")
    } else {
      NULL
    }
  })
  
  
  #TODO it would be nice to include a search box if the number of genes selected exceeds 3 pages
  output$gene_expression_plot<-renderUI({
    if(login$login){
      genes_data<-reactive({
        if(is.null(tissues_samples())){
          #TODO add conditional panel with a value box
          expression_data<-NULL
        } else {
          if(input$gene_tpm_reads=="TPM"){
            table="gene_tpm"
          } else {
            table="gene_reads"
          }
          expression_data<-get_expression(conn=conn, 
                                          gene_id=genes()$geneid, 
                                          tissues_samples=tissues_samples(),
                                          table=table, extra_columns="gene_name"
          )
        }
        return(expression_data)
      })
      
      plot_page<-reactiveValues(current_page=1, max_page=ceiling(length(genes()$geneid)/5))
      
      observeEvent(input$prev_page, {
        if(plot_page$current_page>1){
          plot_page$current_page<-plot_page$current_page-1
        }
      })
      
      observeEvent(input$next_page, {
        if(plot_page$current_page<plot_page$max_page){
          plot_page$current_page<-plot_page$current_page+1
        }
      })
      
      output$page_number<-renderText({
        paste("Page", plot_page$current_page, "of", plot_page$max_page)
      })
      
      plot_data<-reactive({
        if(is.null(genes_data())){
          NULL
        } else {
          page_genes<-unique(genes_data()$gene_id)
          page_genes<-page_genes[((5*(plot_page$current_page-1))+1):(5*plot_page$current_page)]
          plot_data_df<-genes_data() %>% 
            filter(gene_id %in% page_genes)
          return(plot_data_df)
        }
      })
      
      output$gene_boxplot<-renderPlotly({
        if(is.null(plot_data())){
          NULL
        } else {
          plot_ly(data = plot_data(), x=~gene_name, y=~value, 
                  color=~tissue, type = "box", source = "gene_boxplot")%>%
            layout(boxmode = "group", xaxis=list(title="Gene Expression"), 
                   yaxis=list(title=input$gene_tpm_reads)
            )
        }
      })
      
      output$page_navigation<-renderUI({
        if(is.null(genes_data())){
          NULL
        } else {
          tagList(
            br(),
            splitLayout(cellWidths = c("20%", "60%", "20%"),
                        actionLink(inputId = session$ns("prev_page"), 
                                   label=NULL, icon=icon("chevron-left")), 
                        textOutput(outputId = session$ns("page_number")), 
                        actionLink(inputId = session$ns("next_page"), 
                                   label=NULL, icon=icon("chevron-right"))
            ))
        }
      })
      
      tagList(
        plotlyOutput(session$ns("gene_boxplot")),
        column(width = 3, offset = 9, 
               uiOutput(session$ns("page_navigation"))
        )
      )
    } else {
      NULL
    }
  })
}