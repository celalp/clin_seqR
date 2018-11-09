########## FUNCTIONS ###########

get_samples_subjects<-function(datadir, tissues){
  samples_subjects<-list()
  for(tissue in tissues){
    name<-stripWhitespace(tissue)
    name<-gsub(" ", "_", name)
    db_loc<-paste0(datadir, "databases/", name, "_gtex.db")
    con<-dbConnect(drv = SQLite(), dbname=db_loc)
    samples<-dbGetQuery(con, "select * from samples")
    samples$tissue<-tissue
    subjects<-dbGetQuery(con, "select * from subjects")
    subjects$tissue<-tissue
    samples_subjects[["samples"]][[tissue]]<-samples
    samples_subjects[["subjects"]][[tissue]]<-subjects
  }
  samples_subjects$samples<-do.call("rbind", samples_subjects$samples)
  subjid<-strsplit(samples_subjects$samples$SAMPID, split = "-")
  samples_subjects$samples$subject<-paste0("GTEX-", unlist(lapply(subjid, "[", 2)))
  
  samples_subjects$subjects<-do.call("rbind", samples_subjects$subjects)
  samples_subjects$subjects$SEX<-ifelse(samples_subjects$subjects$SEX==1, "Male", "Female")
  samples_subjects$subjects$DTHHRDY<-as.character(samples_subjects$subjects$DTHHRDY)
  samples_subjects$subjects$DTHHRDY<-gsub("0", "Ventilator Case", samples_subjects$subjects$DTHHRDY)
  samples_subjects$subjects$DTHHRDY<-gsub("1", "Fast and Violent", samples_subjects$subjects$DTHHRDY)
  samples_subjects$subjects$DTHHRDY<-gsub("2", "Fast", samples_subjects$subjects$DTHHRDY)
  samples_subjects$subjects$DTHHRDY<-gsub("3", "Intermediate", samples_subjects$subjects$DTHHRDY)
  samples_subjects$subjects$DTHHRDY<-gsub("4", "Slow", samples_subjects$subjects$DTHHRDY)
  samples_subjects$subjects$DTHHRDY[which(is.na(samples_subjects$subjects$DTHHRDY))]<-"Unavailable"
  
  return(samples_subjects)
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
             )#, 
             #tabPanel(title = "Sample Filers", 
            #          tagList(
            #            uiOutput(ns("sample_filters"))
            #          )
            # )
           ) 
           #actionButton(ns("apply_filters"), "Apply Selected Filters")
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
filter_modal_server<-function(input, output, session, tissues, datadir){
  
  output$sample_filters<-renderUI({
    tagList(
      tags$h4("Sample isolation"),
      sliderTextInput(inputId = "autolys",
                      label = "Minimum Autolysis Level:", 
                      choices = c("None", "Mild", "Moderate", "Severe")
      ), 
      sliderTextInput(inputId = "rin",
                      label = "Minimum RNA Integrity (RIN):", 
                      choices = c(1:10),
                      grid = TRUE
      ), 
      sliderTextInput(inputId = "isch",
                      label = "Minimum Ischemic Time:", 
                      choices = c(1:10), #not correct
                      grid = TRUE
      ), 
      sliderTextInput(inputId = "fixative",
                      label = "Minimum time spent on Fixative:", 
                      choices = c(1:10), #not correct
                      grid = TRUE
      ),
      tags$h4("Data Processing"), 
      sliderTextInput(inputId = "map_rate",
                      label = "Minimum mapping rate (%):", 
                      choices = c(1:100),
                      grid = TRUE
      ),
      sliderTextInput(inputId = "exon_rate",
                      label = "Minimum Exon Mapping rate (%):", 
                      choices = c(1:100),
                      grid = TRUE
      ),
      sliderTextInput(inputId = "gene_det",
                      label = "Minimum Genes Detected:", 
                      choices = c(1:10), #not correct
                      grid = TRUE
      ),
      sliderTextInput(inputId = "rin",
                      label = "Minimum Transcripts Detected:", 
                      choices = c(1:10), #not correct
                      grid = TRUE
      )
    )
  })
  
  
  samp_sub<-reactive({
    samp_sub<-get_samples_subjects(datadir = datadir, tissues = tissues())
    #probably needs to be refactored
    
    # these are all subject filters
    if(!is.null(input$sex_filter)){
      samp_sub$subjects<- samp_sub$subjects %>%
        filter(SEX %in% c(input$sex_filter))
    }
    if(!is.null(input$age_filter)){
      samp_sub$subjects<- samp_sub$subjects %>%
        filter(AGE %in% c(input$age_filter))
    }
    if(!is.null(input$death_filter)){
      samp_sub$subjects<- samp_sub$subjects %>%
        filter(DTHHRDY %in% c(input$death_filter))
    }
    if(input$all_or_common == "All"){
      NULL
    } else {
      common_subjects<-split(samp_sub$subjects, samp_sub$subjects$tissue)
      common_subjects <- common_subjects %>%
        #######################################################################
      purrr::reduce(inner_join, by="SUBJID") #this is not working as expected
      #######################################################################
      subj<-common_subjects$SUBJID
      samp_sub$subjects<- samp_sub$subjects %>%
        filter(SUBJID %in% subj)
      # sample filters go here
    }
    
    
    #finally filter remaining samples by subject id
    
    samp_sub$samples<-samp_sub$samples %>%
      filter (subject %in% c(samp_sub$subjects$SUBJID))
    return(samp_sub)
  })
  
  
  
  output$sex_dist<-renderPlotly({
    for_bar<-plyr::count(samp_sub()$subjects[,c(3,6)])
    plot_ly(data = for_bar, x=~tissue, y=~freq, color=~SEX, type='bar') %>%
      layout(yaxis = list(title = 'Count'), barmode = 'stack')
  })
  
  output$age_dist<-renderPlotly({
    for_bar<-plyr::count(samp_sub()$subjects[,c(4,6)])
    for_bar<-for_bar[base::order(for_bar$AGE),]
    plot_ly(data = for_bar, x=~tissue, y=~freq, color=~AGE, type='bar') %>%
      layout(yaxis = list(title = 'Count'), barmode = 'stack')
  })
  
  return(reactive({
    tissues<-unique(samp_sub()$samples$tissue)
    for_server<-split(samp_sub()$samples[,c(2,65)], f = samp_sub()$samples$tissue)
    names(for_server)<-tissues
    return(for_server)
  }))
}





