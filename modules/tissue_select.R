tissue_select_ui<-function(id){
  ns<-NS(id)
  tagList(
    column(width = 4,
           tagList(
             uiOutput(ns("selected_tissues_ui")),
             uiOutput(ns("filters"))
           )
    ),
    column(width=8, 
           tagList(
             uiOutput(ns("plots"))
           )
    )
  )
}

tissue_select<-function(input, output, session, conn, login){
  output$selected_tissues_ui<-renderUI({
    if(login$login){
      query<-"select distinct(smtsd) from samples.samples;"
      query<-sqlInterpolate(conn, query)
      choices<-dbGetQuery(conn, query)[,1]
      choices<-unlist(gsub("_", " ", choices))
      choices<-choices[order(choices)]
      tagList(
        pickerInput(session$ns("selected_tissues"), label = "Select Tissues",
                    choices = choices, 
                    options=list(`max-options`=10, 
                                 `live-search`=T, size=10), 
                    multiple = TRUE)
      )
    } else {
      NULL
    }
  })
  
  output$filters<-renderUI({
    if(login$login){
      tagList(
        tags$h2("Subject Filters"),
        tagList(
          radioGroupButtons(session$ns("all_or_common"), label = "Use all or overlapping samples", 
                            choices = c("All", "Overlapping"), direction = "vertical", 
                            selected = "All"), 
          checkboxGroupInput(inputId = session$ns("sex_filter"), label = "Sex", 
                             choices = c("Male", "Female")), 
          checkboxGroupInput(inputId = session$ns("age_filter"), label = "Age group", 
                             choices = c("20-29", "30-39", "40-49", 
                                         "50-59", "60-69", "70-79")),
          checkboxGroupInput(inputId = session$ns("death_filter"), label = "Cause of death", 
                             choices = c("Ventilator Case", "Fast and Violent", "Fast", 
                                         "Intermediate", "Slow", "Unavailable"))
        )
      )
      # actionButton(session$ns("apply_filter"), label = "Apply Filters")
    } else {
      NULL
    }
  })
  
  samp_sub<-reactive({
    samp_sub_df<-get_samples_subjects(conn = conn, tissues = input$selected_tissues)
    if(!is.null(samp_sub_df)){
      if(!is.null(input$sex_filter)){
        samp_sub_df<- samp_sub_df %>%
          filter(sex %in% c(input$sex_filter))
      }
      if(!is.null(input$age_filter)){
        samp_sub_df<- samp_sub_df %>%
          filter(age %in% c(input$age_filter))
      }
      if(!is.null(input$death_filter)){
        samp_sub_df<- samp_sub_df %>%
          filter(dthhrdy %in% c(input$death_filter))
      }
      if(input$all_or_common == "All"){
        NULL
      } else {
        common_subj<-get_common(samp_sub_df)
        samp_sub_df<- samp_sub_df %>%
          filter(subjid %in% common_subj)
      }
      
      #TODO add sample filters here
      return(samp_sub_df)
    } else {
      NULL
    }
  })
  
  output$plots<-renderUI({
    if(login$login){
      output$sex_dist<-renderPlotly({
        if(is.null(samp_sub())){
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
      
      output$tissue_alert<-renderUI({
        alert_box(alert = "Please select tissues to see population distribution", color = "maroon", 
                  condition = is.null(samp_sub()))
      })
      
      tagList(
        uiOutput(session$ns("tissue_alert")),
        tags$h4("Sex Distribution"),
        plotlyOutput(session$ns("sex_dist")), 
        br(),
        tags$h4("Age Distribution"),
        plotlyOutput(session$ns("age_dist"))
      )
    } else {
      NULL
    }
  })
  
  return(reactive({
    if(!is.null(samp_sub())){
      tissues_samples<-get_tissues_samples(samp_sub())
      return(tissues_samples)
    } else {
      NULL
    }
  }))
}

