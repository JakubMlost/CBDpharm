library(shiny)
library(lobstr)
library(tibble)
library(VennDiagram)
library(RJSONIO)
library(STRINGdb)
library(feather)
library(grid)
library(futile.logger)
library(influential)
library(igraph)

#some code from stackoverflow, which allows the app to deploy 
options(repos = BiocManager::repositories())
getOption("repos")

#Download string DB
string_db <- STRINGdb$new( version="11", species=9606, score_threshold=400, input_directory="")
#Download protein names for string IDs
terms<-string_db$get_aliases()
#Load list of diseases on OpenTargets
diseasePath <- "20.11_disease_list.feather"
dt<-read_feather(diseasePath)
dt<-as.data.frame(dt)
diseases<-dt[,2]

#Load CBD data
CBDPath <- "net_below_2mM.csv"
CBDallPath <- "net.csv"
CBDall<-read.table(CBDallPath, sep=";", header=TRUE)
CBD<-read.table(CBDPath, sep=";", header=TRUE)
CBDtargets<-unique(CBD$ensembl, na.rm=TRUE)
CBDtarget <- CBDtargets[!is.na(CBDtargets)]
CBDtargets<-c(CBDtarget)
myHits <- read.table("CBDstringID.csv", header=TRUE)
myHits<-c(myHits$stringId)
#UI
ui <- fluidPage(
                  navlistPanel("CBDpharm",
                               tabPanel(
                                 selectizeInput('foo', label = "What is your interest?", choices = c("Please enter name",diseases),multiple = FALSE),
                                 helpText(tags$strong("Welcome to CBDpharm!"),br(),
                                        "Please enter the name of the disease you are interested in on the left and give us ~15s to compute."),br(),
                                        "P.S. First, you need to delete the current text in the select input box!"),
                               tabPanel("CBD pharmacology",
                                        dataTableOutput("CBD")),
                               tabPanel("Common Targets",
                                        uiOutput(outputId = "image")),
                               tabPanel("Enrichment Analysis",
                                        dataTableOutput("enrich")),         
                               tabPanel("Influential Node",
                                        dataTableOutput("test")),
                               tabPanel("Influential graph",
                                        wellPanel(fluidRow(
                                          column(3,
                                                   sliderInput(
                                                     "min", label = "min.node.size:",
                                                     min = 1, value = 20, max = 50
                                                   )),
                                          column(3,
                                                 sliderInput(
                                                   "max", label = "max.node.size:",
                                                   min = 10, value = 30, max = 100
                                                 )),
                                          column(2,
                                                 sliderInput(
                                                   "distpower", label = "dist.power:",
                                                   min = 0.1, value = 0.5, max = 1
                                                 )),
                                          column(2,
                                                 sliderInput(
                                                   "labelcex", label = "label.cex:",
                                                   min = 0.1, value = 0.15, max = 1
                                                 )),
                                          column(1,
                                                 selectizeInput("color", label = "Color", choices = c("viridis","magma","inferno","plasma","cividis")),
                                                 ),
                                          column(1,
                                                 selectizeInput("layout", label = "Layout", choices = c("sphere","kk", "star", "components", "circle", "automatic", "grid", "random", "dh", "drl", "fr", "gem", "tree","graphopt", "lgl", "mds"), selected="kk"),
                                              )),
                                        fluidRow(plotOutput("influential",   width = "100%",
                                                            height = "800px",)),
                                        )),
                               tabPanel("About", helpText(tags$p("CBD pharm is a web application predestined for preliminary analysis of CBD's molecular targets in disease of interest (DOI)."),
                                                                  tags$p("The app uses self-curated database of CBD's molecular targets (Described in Mlost et al., 2020 in ", tags$cite("International Journal of Molecular Science"),")",
                                                                  tags$p("Venn analysis to detect common targets between CBD and DOI"), 
                                                                  tags$p("string-db to visualise common targets network and analyse their protein-protein interaction with targets within DOI netwok"),
                                                                  tags$p("KEGG database for functional enrichment analysis of common targets, which reveal the biological processes involved in CBD's mechanism of action"),
                                                                  tags$p("The final outcome of the analysis is the list of CBD's molecular targets that has the highest number of connections 
                                                                         with DOI targets (degree of the node) combined with the CBD's affinity/IC50/EC50 value [nM] and interaction type with the given target."),
                                                                  br(),br(),"To use the app, first enter the name of the disease you are interested in.",
                                                                                                    "After ~15s computations should be done and you can view your results."))),
                               tabPanel("References", helpText(tags$ol(tags$li(
                                 "Szklarczyk, Damian, et al.", tags$cite("STRING v11: protein – protein association networks with increased coverage, supporting functional discovery in genome - wide experimental datasets."), "Nucleic acids research 47.D1 (2019): D607 - D613."
                                ),
                                tags$li("Kanehisa, Minoru, and Susumu Goto.", tags$cite("KEGG: kyoto encyclopedia of genes and genomes."),
                                "Nucleic acids research 28.1 (2000): 27-30."),
                                tags$li("Mlost, Jakub, Marta Bryk, and Katarzyna Starowicz.", tags$cite("Cannabidiol for Pain Treatment: Focus on Pharmacology and Mechanism of Action."),
                                "International journal of molecular sciences 21.22 (2020): 8870.")))),  
                               fluid = FALSE, widths = c(2, 10))
)
                  
#SERVER
server <- function(input, output, session) {
  
  output$CBD<-renderDataTable(
    CBDall,   options = list(
      pageLength = 15)
  )
  
  
  #create link to diseases
  index<-reactive({
    req(input$foo != "")
    index<-which(dt[,2] == input$foo)
    URL<-paste("https://platform-api.opentargets.io/v3/platform/public/association/filter?disease=", dt[index,1],"&size=10000&fields=target.id", sep="")
  })
  
  #create table with targets of interest
  diseaseTargets<- reactive({
    foodMarketsRaw<-RJSONIO::fromJSON(index())
    foodMarkets<-foodMarketsRaw[['data']]
    foodMarkets[[1]][[1]][[1]]
    targetID<-sapply(foodMarkets, function(x) x[[1]][[1]])
  })
  
  #Venn analysis
  VennTable<-reactive({
    x<-list(diseaseTargets(),CBDtargets)
    y<-calculate.overlap(x)
    VennTable<-y$a3
  })
  
  
  #Render enrichment results
  output$enrich<-renderDataTable({
    validate(
      need(VennTable(), 'There are no common targets!')
    )
    id <- showNotification("Please wait while the data is loading...", duration = NULL, closeButton = FALSE)
    on.exit(removeNotification(id), add = TRUE)
    enrichment <- string_db$get_enrichment(VennTable(), category = "KEGG")
    x<-head(enrichment, n=20)
    y<-subset(x, select = c(description, p_value, fdr, preferredNames))},list(
      pageLength = 10))
    #plot Venn 
  output$image <- renderUI({
    validate(
      need(VennTable(), 'There are no common targets!')
    )
    string<-paste(VennTable(),collapse = "%0d")
    URL = paste("https://string-db.org/api/image/network?identifiers=",string,"&species=9606",sep="")
    tags$img(src = URL)
  })    #Plot string
    interactions<-reactive({
      fmNames<-diseaseTargets()
      source<-c()
      target<-c()
      score<-c()
      min<-1
      max<-100
      while (min<length(fmNames)){
        string<-paste(fmNames[min:max],collapse = "%0d")
        URL = paste("https://string-db.org/api/tsv/network?identifiers=",string,"&species=9606",sep="")
        table<-read.table(URL, sep="\t", header=TRUE)
        tempSource<-c(table$stringId_A)
        source<-c(source, tempSource)
        tempTarget<-c(table$stringId_B)
        target<-c(target,tempTarget)
        tempScore<-c(table$score)
        score<-c(score, tempScore)
        min<-min+100
        max<-max+100
      }
      interactions<-cbind.data.frame(source,target,score)
  }) 
      output$test<-renderDataTable({
        validate(
          need(VennTable(), 'There are no common targets!')
        )
        id <- showNotification("Please wait while the data is loading...", duration = NULL, closeButton = FALSE)
        on.exit(removeNotification(id), add = TRUE)
      #get interactions of interest
      myInteractions=c()
      for (i in myHits){
        tempTable <- subset(interactions(), source==i | target==i)
        Freq<-c(length(tempTable$source))
        Var1<-c(i)
        tempTable<-data.frame(Var1,Freq)
        myInteractions <- rbind(myInteractions, tempTable) 
      } 
      newdata <- myInteractions[order(-myInteractions$Freq),] 
      
      #rename string IDs to protein names
      x=1
      activity<-c()
      interaction<-c()
      for(i in newdata$Var1){
        index<-which(terms$STRING_id==i)
        index<-index[1]
        newdata$Var1[x]=terms$alias[index]
        index<-which(CBD$gene==newdata$Var1[x])
        if(is.na(CBD$Affinity.Ki..nM.[index][1])){
          activity<-c(activity, min(CBD$EC.IC50..nM.[index],na.rm=T))}
        else{activity<-c(activity,min(CBD$Affinity.Ki..nM.[index],na.rm=T))}
        if (is.na(CBD$Interaction.type[index[1]])){interaction<-c(interaction, CBD$Interaction.type[index[2]])}
        else {interaction<-c(interaction, CBD$Interaction.type[index[1]])}
        x=x+1
      }
      activity<-as.numeric(activity)
      newdata<-cbind(newdata,activity,interaction)
      colnames(newdata)<-c("Target","Vertex degree","Activity[nM]","interaction")
      newdata
    },list(
      pageLength = 15)
)
      My_graph<-reactive({
        id <- showNotification("Please wait while the data is loading...", duration = NULL, closeButton = FALSE)
        on.exit(removeNotification(id), add = TRUE)
        #get interactions of interest
        MyData=c()
        for (i in myHits){
          tempTable <- subset(interactions(), source==i | target==i)
          MyData <- rbind(MyData, tempTable) 
        } 
        x=1
        for(i in MyData[,1]){
          index<-which(terms$STRING_id==i)
          index<-index[1]
          MyData[x,1]<-terms$alias[index]
          x=x+1
        }   
        
        # Preparing the data
        x=1
        for(i in MyData[,2]){
          index<-which(terms$STRING_id==i)
          index<-index[1]
          MyData[x,2]<-terms$alias[index]
          x=x+1
        } 
        My_graph <- graph_from_data_frame(d=MyData)
      })
  
    output$influential<-renderPlot({
      validate(
        need(VennTable(), 'There are no common targets!')
      )
      # Extracting the vertices
      GraphVertices <- V(My_graph())   
      # Calculating degree centrality
      My_graph_degree <- degree(My_graph(), v = GraphVertices, normalized = FALSE)
      # Visualizing the graph based on IVI values
      My_graph_IVI_Vis <- cent_network.vis(graph = My_graph(),
                                           node.color=input$color,
                                           dist.power=input$distpower,
                                           label.cex=input$labelcex,
                                           cent.metric = My_graph_degree,
                                           directed = FALSE,
                                           plot.title = "Node centrality degree",
                                           legend.title = "value",
                                           layout=input$layout,
                                           stroke.size=1,
                                           node.size.min=input$min,
                                           node.size.max=input$max,
                                           show.bottom.border=FALSE,
                                           show.left.border=FALSE
                                           
                                           
                                           
      )
      My_graph_IVI_Vis
      })
} 
  shinyApp(ui = ui, server = server)
