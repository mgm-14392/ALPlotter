library(shiny)

shinyUI(fluidPage(
  titlePanel("Activity Landscape Plotter V.1"),
  tabsetPanel(
    tabPanel("About",
     sidebarPanel(h2("DIFACQUIM"),
                  br(),
                  p("Welcome to this ")
                  
                  ),
     mainPanel(
       h2("What are they for?")
     )
    ),
    tabPanel("Instructions",
             sidebarPanel(
               
             ),
    mainPanel(
      
    )
      
    ),
    tabPanel("SAS map",
  sidebarLayout(
    sidebarPanel(
      fileInput("file1", "Upload the csv file",   
                accept = c("text/csv",
        "text/comma-separated-values,text/plain",
        ".csv")),
      selectInput("FP",
                  "Select a fingerprint",
                  c("ECFP","Pubchem","MACCS"
                  )),
      radioButtons("dep",
                   "Select the diameter for topological fingerprints",
                   choices=c("4"=4, "6"=6)),
      selectInput("act","Select the column with the activity you want to use", 
                  choices = c("Activity 1"= 2, "Activity 2"= 3, "Activity 3"= 4, "Activity 4"= 5)),
      selectInput("col","Select the color for the data points", 
                  choices = c("SALI", "Max.Activity")),
      numericInput("ablineX", "Choose the X axis threshold",
                   min = 0, max = 1, value =0.5, step = 0.5),
      numericInput("ablineY", "Choose the Y axis threshold",
                   min = 0, max = 5, value = 2, step = 0.5),
      downloadButton('downloadData', 'Download SASmap data')
    
    ),
    mainPanel(
      plotOutput("sasplot", width = 800, height = 600,  brush = brushOpts(id = "plot_brush")),
      column(width = 12,
             h4("You selected"),
             verbatimTextOutput("brush_info")
      ),
      downloadButton("down","Download image"),
              downloadButton('ACdata', 'Activity Cliffs'),
              downloadButton('C1data', 'Non Descriptive'),
              downloadButton('SSARdata', 'Smooth SAR'),
              downloadButton('RHata', 'Similarity Cliffs'),
              downloadButton('UX', 'On the X axis threshold'),
              downloadButton('UY', 'On the Y axis threshold'))
  )
    ),
  tabPanel("Density SAS map",
           sidebarLayout(
             sidebarPanel(
               fileInput("file3", "Upload the csv file",   
                         accept = c("text/csv",
                                    "text/comma-separated-values,text/plain",
                                    ".csv")),
               selectInput("FP3",
                           "Select a fingerprint",
                           c("ECFP","Pubchem","MACCS"
                           )),
               radioButtons("dep3",
                            "Select the diameter for topological fingerprints",
                            choices=c("4"=4, "6"=6)),
               selectInput("act4","Select the column with the activity you want to use", 
                           choices = c("Activity 1"= 2, "Activity 2"= 3, "Activity 3"= 4, "Activity 4"= 5)),
               numericInput("ablineX2", "Choose the X axis threshold",
                            min = 0, max = 1, value =0.5, step = 0.5),
               numericInput("ablineY2", "Choose the Y axis threshold",
                            min = 0, max = 5, value = 2, step = 0.5),
               downloadButton('downloadData2', 'Download SASmap data')
               
             ),
             mainPanel(
               plotOutput("sasplot2", width = 800, height = 600,  brush = brushOpts(id = "plot_brush2")),
               column(width = 12,
                      h4("You selected"),
                      verbatimTextOutput("brush_info2")
               ),
               downloadButton("down2","Download image"),
               downloadButton('ACdata2', 'Activity Cliffs'),
               downloadButton('C1data2', 'Non Descriptive'),
               downloadButton('SSARdata2', 'Smooth SAR'),
               downloadButton('RHata2', 'Similarity Cliffs'),
               downloadButton('UX2', 'On the X axis threshold'),
               downloadButton('UY2', 'On the Y axis threshold'))
           )
  ),
  
  tabPanel("DAD map",
           sidebarLayout(
             sidebarPanel(
               fileInput("file2", "Upload the csv file",   
                         accept = c("text/csv",
                                    "text/comma-separated-values,text/plain",
                                    ".csv")),
               selectInput("FP2",
                           "Select a fingerprint to compute the similarity and add color to the DAD map",
                           c("ECFP","Pubchem","MACCS"
                           )),
               radioButtons("dep2",
                            "Select the diameter for topological fingerprints",
                            choices=c("4"=4, "6"=6)),
               selectInput("act2", "Select the column with the activity you want to use",
                           choices = c("Activity 1"=2, "Activity 2"=3,"Activity 3"=4, "Activity 4"=5)),
               selectInput("act3", "Select the column with the activity you want to use",
                           choices=c("Activity 2"=3,"Activity 1"=2, "Activity 3"=4, "Activity 4"=5)),
               selectInput("colp","Select the color for the data points", 
                           choices = c("Similarity", "Selectivity")),
               numericInput("minlinex", "Choose the lowest X axis threshold",
                            min= -100, max= 100, value= -2, step= 0.5),
               numericInput("maxlinx", "Choose the highest X axis theshold",
                            min = -100, max = 100, value = 2, step = 0.5),
               numericInput("minliney", "Choose the lowest Y axis threshold",
                            min = -100, max = 100, value = -2, step = 0.5),
               numericInput("maxliney", "Choose the highest Y axis threshold",
                            min = -100, max = 100, value = 2, step = 0.5),
               downloadButton('downloadDADData', 'Download DADmap data')
             ),
             mainPanel(plotOutput("DADmap",
                                  width = 800, height = 600, brush = brushOpts(id = "plot_brush1")),
                       column(width = 12,
                              h4("You selected"),
                              verbatimTextOutput("brush_info1")
                                  ),
                       downloadButton("down3","Download image"),
                       downloadButton('ZON1S', 'Z1u'),
                       downloadButton('ZON2I', 'Z2d'),
                       downloadButton('ZON2S', 'Z2u'),
                       downloadButton('ZON1I', 'Z1d'),
                       downloadButton('ZON3S', 'Z3u'),
                       downloadButton('ZON3I', 'Z3d'),
                       downloadButton('ZON4L', 'Z4l'),
                       downloadButton('ZON4R', 'Z4r'),
                       downloadButton('ZON5', 'Z5'))
           )
           ),
  
  tabPanel("Density DAD map",
           sidebarLayout(
             sidebarPanel(
               fileInput("file4", "Upload the csv file",   
                         accept = c("text/csv",
                                    "text/comma-separated-values,text/plain",
                                    ".csv")),
               selectInput("act5", "Select the column with the activity you want to use",
                           choices = c("Activity 1"=2, "Activity 2"=3,"Activity 3"=4, "Activity 4"=5)),
               selectInput("act6", "Select the column with the activity you want to use",
                           choices=c("Activity 2"=3,"Activity 1"=2, "Activity 3"=4, "Activity 4"=5)),
               numericInput("mlx", "Choose the lowest X axis threshold",
                            min= -100, max= 100, value= -2, step= 0.5),
               numericInput("maxlx", "Choose the highest X axis theshold",
                            min = -100, max = 100, value = 2, step = 0.5),
               numericInput("mley", "Choose the lowest Y axis threshold",
                            min = -100, max = 100, value = -2, step = 0.5),
               numericInput("maxly", "Choose the highest Y axis threshold",
                            min = -100, max = 100, value = 2, step = 0.5),
               downloadButton('downloadDADData3', 'Download DADmap data')
             ),
             mainPanel(plotOutput("DADmap2",
                                  width = 800, height = 600, brush = brushOpts(id = "plot_brush3")),
                       column(width = 12,
                              h4("You selected"),
                              verbatimTextOutput("brush_info3")
                       ),
                       downloadButton("down4","Download image"),
                       downloadButton('ZON1S2', 'Z1u'),
                       downloadButton('ZON2I2', 'Z2d'),
                       downloadButton('ZON2S2', 'Z2u'),
                       downloadButton('ZON1I2', 'Z1d'),
                       downloadButton('ZON3S2', 'Z3u'),
                       downloadButton('ZON3I2', 'Z3d'),
                       downloadButton('ZON4L2', 'Z4l'),
                       downloadButton('ZON4R2', 'Z4r'),
                       downloadButton('ZON52', 'Z5'))
           )
  )
  
  )
)
)