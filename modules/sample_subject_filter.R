########## FUNCTIONS ###########

get_samples_subjects<-function(conn, tissues){
  if(!is.null(tissues)){
    samples_subjects<-list()
    for(tissue in tissues){
      name<-gsub(" ", "_", tissue)
      sample_table<-paste(name, "samples", sep=".")
      subject_table<-paste(name, "subjects", sep=".")
      sample_query<-sqlInterpolate(conn, paste("select * from", sample_table))
      subject_query<-sqlInterpolate(conn, paste("select * from", subject_table))
      samples<-dbGetQuery(conn, sample_query)
      samples$tissue<-tissue
      subjects<-dbGetQuery(conn, subject_query)
      subjects$tissue<-tissue
      samples_subjects[["samples"]][[tissue]]<-samples
      samples_subjects[["subjects"]][[tissue]]<-subjects
    }
    samples_subjects$samples<-do.call("rbind", samples_subjects$samples)
    samples_subjects$subjects<-do.call("rbind", samples_subjects$subjects)
    subjid<-strsplit(samples_subjects$samples$sampid, split = "-")
    samples_subjects$samples$subject<-paste0("GTEX-", unlist(lapply(subjid, "[", 2)))
    samples_subjects$subjects$sex<-ifelse(samples_subjects$subjects$sex==1, "Male", "Female")
    samples_subjects$subjects$dthhrdy<-as.character(samples_subjects$subjects$dthhrdy)
    samples_subjects$subjects$dthhrdy<-gsub("0", "Ventilator Case", samples_subjects$subjects$dthhrdy)
    samples_subjects$subjects$dthhrdy<-gsub("1", "Fast and Violent", samples_subjects$subjects$dthhrdy)
    samples_subjects$subjects$dthhrdy<-gsub("2", "Fast", samples_subjects$subjects$dthhrdy)
    samples_subjects$subjects$dthhrdy<-gsub("3", "Intermediate", samples_subjects$subjects$dthhrdy)
    samples_subjects$subjects$dthhrdy<-gsub("4", "Slow", samples_subjects$subjects$dthhrdy)
    samples_subjects$subjects$dthhrdy[which(is.na(samples_subjects$subjects$dthhrdy))]<-"Unavailable"
    return(samples_subjects)
  } else {
    NULL
  }
}



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
             ), 
             actionButton(ns("apply_filters"), "Apply Selected Filters")
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
        samp_sub$subjects<- samp_sub$subjects %>%
          filter(sex %in% c(input$sex_filter))
      }
      if(!is.null(input$age_filter)){
        samp_sub$subjects<- samp_sub$subjects %>%
          filter(age %in% c(input$age_filter))
      }
      if(!is.null(input$death_filter)){
        samp_sub$subjects<- samp_sub$subjects %>%
          filter(dthhrdy %in% c(input$death_filter))
      }
      if(input$all_or_common == "All"){
        NULL
      } else {
        common_subjects<-split(samp_sub$subjects, samp_sub$subjects$tissue)
        common_subjects <- common_subjects %>%
          #######################################################################
        purrr::reduce(inner_join, by="subjid") #this is not working as expected
        #######################################################################
        subj<-common_subjects$subjid
        samp_sub$subjects<- samp_sub$subjects %>%
          filter(subjid %in% subj)
      }
      #this has the sample filters applied
      samp_sub$samples<-samp_sub$samples %>%
        filter (subject %in% c(samp_sub$subjects$subjid))
      #TODO apply sample filters
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
      for_bar<-plyr::count(samp_sub()$subjects[,c(3,6)])
      plot_ly(data = for_bar, x=~tissue, y=~freq, color=~sex, type='bar') %>%
        layout(yaxis = list(title = 'Count'), barmode = 'stack')
    }
  })
  
  output$age_dist<-renderPlotly({
    if(is.null(samp_sub())){
      NULL
    } else {
      for_bar<-plyr::count(samp_sub()$subjects[,c(4,6)])
      for_bar<-for_bar[base::order(for_bar$age),]
      plot_ly(data = for_bar, x=~tissue, y=~freq, color=~age, type='bar') %>%
        layout(yaxis = list(title = 'Count'), barmode = 'stack')
    }
  })
  
  return(reactive({
    if(!is.null(samp_sub)){
      tissues<-unique(samp_sub()$samples$tissue)
      for_server<-split(samp_sub()$samples[,c(2,65)], f = samp_sub()$samples$tissue)
      names(for_server)<-tissues
      return(for_server)
    }
  }))
}





