
########## MODULE ###########

filter_modal_ui<-function(id){
  ns<-NS(id)
  fluidRow(
    column(width = 4,
           tabsetPanel(
             tabPanel(title = "Subject Filters",
                      tagList(
                        radioGroupButtons(ns("all_or_common"), label = "Use all or overlapping samples", 
                                          choices = c("All", "Overlapping"), direction = "vertical", 
                                          selected = "All"), 
                        checkboxGroupInput(inputId = ns("sex_filter"), label = "Sex", 
                                           choices = c("Male", "Female")), 
                        checkboxGroupInput(inputId = ns("age_filter"), label = "Age group", 
                                           choices = c("20-29", "30-39", "40-49", 
                                                       "50-59", "60-69", "70-79")),
                        checkboxGroupInput(inputId = ns("death_filter"), label = "Cause of death", 
                                           choices = c("Ventilator Case", "Fast and Violent", "Fast", 
                                                       "Intermediate", "Slow", "Unavailable"))
                      )
             ),
             tabPanel(title="Sample Filters",
                      tagList(
                        tags$h4("Sample isolation"),
                        sliderTextInput(inputId = ns("autolys"),
                                        label = "Minimum Autolysis Level:", 
                                        choices = c("None", "Mild", "Moderate", "Severe")
                        ), 
                        sliderTextInput(inputId = ns("rin"),
                                        label = "Minimum RNA Integrity (RIN):", 
                                        choices = c(1:10),
                                        grid = TRUE
                        ),
                        tags$h4("Data Processing"), 
                        sliderTextInput(inputId = ns("map_rate"),
                                        label = "Minimum mapping rate (%):", 
                                        choices = c(1:100),
                                        grid = TRUE
                        )
                      )
             )
           )
    ),
    column(width = 8, 
           tagList(
             tags$h4("Sex Distribution"),
             plotlyOutput(ns("sex_dist")), 
             br(),
             tags$h4("Age Distribution"),
             plotlyOutput(ns("age_dist"))
           )
    ) 
  )
}


# tissues will come from input$selected_tissues
filter_modal_server<-function(input, output, session, tissues, conn){
  
  samp_sub<-reactive({
    samp_sub<-get_samples_subjects(conn = conn, tissues = tissues())
    if(!is.null(samp_sub)){
      if(!is.null(input$sex_filter)){
        samp_sub<- samp_sub %>%
          filter(sex %in% c(input$sex_filter))
      }
      if(!is.null(input$age_filter)){
        samp_sub<- samp_sub %>%
          filter(age %in% c(input$age_filter))
      }
      if(!is.null(input$death_filter)){
        samp_sub<- samp_sub %>%
          filter(dthhrdy %in% c(input$death_filter))
      }
      if(input$all_or_common == "All"){
        NULL
      } else {
        common_subj<-get_common(samp_sub)
        samp_sub<- samp_sub %>%
          filter(subjid %in% common_subj)
      }
      
      #TODO add sample filters here
      
      return(samp_sub)
    } else {
      NULL
    }
  })
  
  
  
  output$sex_dist<-renderPlotly({
    if(is.null(samp_sub())){
      toastr_error("Please select tissues to see demographic and sample distribution") 
      NULL
    } else {
      for_bar<-plyr::count(samp_sub()[,c("sex","smtsd")])
      colnames(for_bar)[2]<-"tissue"
      plot_ly(data = for_bar, x=~tissue, y=~freq, color=~sex, type='bar') %>%
        layout(yaxis = list(title = 'Count'), barmode = 'stack')
    }
  })
  
  output$age_dist<-renderPlotly({
    if(is.null(samp_sub())){
      NULL
    } else {
      for_bar<-plyr::count(samp_sub()[,c("age","smtsd")])
      colnames(for_bar)[2]<-"tissue"
      for_bar<-for_bar[base::order(for_bar$age),]
      plot_ly(data = for_bar, x=~tissue, y=~freq, color=~age, type='bar') %>%
        layout(yaxis = list(title = 'Count'), barmode = 'stack')
    }
  })
  
  return(reactive({
    if(!is.null(samp_sub)){
      for_server<-split(samp_sub()[,c("sampid","smtsd")], 
                        f = as.factor(samp_sub()$smtsd))
      names(for_server)<-tissues()
      return(for_server)
    }
  }))
}





