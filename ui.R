library(shiny)

shinyUI(fluidPage(
  titlePanel("Activity Landscape Plotter V.1"),
  tabsetPanel(
    tabPanel("About",
     sidebarPanel(h2("DIFACQUIM"),
                  br(),
                  p("Welcome to this first online
                    Activity Landscape Plotter. Here you will be able to perform
                    analyses of activity landscape using Structure-Activity Similarity
                    (SAS) maps, Dual Activity-Difference (DAD) maps and other metrics.",style="text-align:justify"),
                  p("If you use this app please cite this reference:"),
                  p(a("Gonzalez-Medina M, Mendez-Lucio O, Medina-Franco JL. Activity Landscape Plotter: A Web-based Application for the Analysis of Structure-Activity Relationships.",
                    em("J Chem Inf Model"),"(2017) 57(3), 397-402.", href="http://pubs.acs.org/doi/abs/10.1021/acs.jcim.6b00776", 
                    target="_blank", style="text-align:justify")),
                  p("This App was developed by", a("Mariana Gonzalez-Medina", 
                                                   href="https://www.researchgate.net/profile/Mariana_Gonzalez-Medina2", target="_blank"),
                  a("Oscar Mendez-Lucio",  href="https://www.researchgate.net/profile/Oscar_Mendez-Lucio", target="_blank"), "and", 
                    a("Jose L. Medina-Franco", href="https://www.researchgate.net/profile/Jose_Medina-Franco/publications", target="_blank"), "members of", style="text-align:justify", 
                  a("DIFACQUIM", href="http://www.difacquim.com/", target="_blank"),
                    "The group is based in the Pharmacy Department of the School of Chemistry in", 
                    a("Universidad Nacional Autonoma de Mexico.", href="http://english.unam.mx/", target="_blank")),
                  img(src="Logo.png", height=150, width=120,
                      style="display: block; margin-left: auto; margin-right: auto"),
                  br(),
                  strong("Last update July, 2017")
                  
                  ),
     mainPanel(
       h2("What are they for?"),
       h4("Structure-Activity Similarity (SAS) maps."),
       p("SAS maps were introuced to find a relationship between structure and activity, 
         based on a systematic pairwise comparison of all the compounds in a data set. You can find more information",
         style="text-align:justify",
         a("here.", href="http://pubs.acs.org/doi/abs/10.1021/ci300362x", target="_blank")),
       p("With this App you can automatically obtain all the information required to generate and analyze SAS maps
         using your own activity data for one or more biological targets. You can customize the SAS maps and
         download the raw data as described in the Instructions."),
       h4("Dual Activity-Difference (DAD) maps."),
       p("DAD maps depict pairwise activity differences for 
         each possible pair of compounds in a data set against two targets. These maps are helpful to differentiate
         when a strucutural modification increases or decreases the activity for one target or the other. You can find more
         information", style="text-align:justify", a("here.", 
                                     href="http://pubs.rsc.org/en/Content/ArticleLanding/2011/MD/c0md00159g#!divAbstract",
                                     target="_blank")),
       p("With this App you can automatically obtain all the information required to plot and analyze DAD maps,
         using the biological activities in your input file.
         You can customize the DAD maps and download the raw data as described in the Instructions.")
       
     )
    ),
    tabPanel("Instructions",
             sidebarPanel(
               p(strong("Must of the errors occur beacause your input file does not have the right format. See the template if you are not sure about the format.", style="text-align:justify"),style="text-align:justify"),
               p(span("1.", style="color:red"),"How to introduce your data:"),
               p("You will need a comma delimited file with these columns:", style="text-align:justify"),
               p("First column: SMILES, do not leave empty rows", style="text-align:justify"),
               p("Second column: first activity, it is okay if you leave empty rows in this column", style="text-align:justify"),
               p("Third column: second activity, same as for the second column. If you do not have more than one 
                 activity leave the rest of the activity columns empty", style="text-align:justify"),
               p("Fourth column: third activity, leave this column empty if you do not have a third activity", style="text-align:justify"),
               p("Fifth column: Same as for the fourth and third columns", style="text-align:justify"),
               p("Sixth column:",strong("this column is for your compounds IDs, make sure the last column in your file has the IDs.",
                                        style="text-align:justify"), style="text-align:justify"),
               p(strong("You have to introduce your activity data as PIC50 =-log10(IC50) or pkI=-log10(kI50).",
                      style="text-align:justify"), style="text-align:justify"),
               br(),
p("You can download the template from the template tab."),
p(span("2.", style="color:red"), "After you upload your data, a plot with the pre-stablished settings on the App will appear.
  This could take a few seconds if you are analyzing many compounds.", style="text-align:justify", strong("Wait until the first
                                                                              plot appears before you start changing the settings.")),
p(span("3.", style="color:red"), "You can choose a fingerpirnt, the column with activity and adjust the thresholds.", style="text-align:justify"),
p(span("4.", style="color:red"), "You can choose how to color the data points (the most active compound in the pair, SALI value, selectivity or similarity).", style="text-align:justify"),
p(span("5.", style="color:red"), "Whenever you change an input,",
  strong("wait until the plot is completely loaded before you download new data.", style="text-align:justify"))
               
             ),
    mainPanel(
      h3("Information about each plot and the metrics used."),
      h4(strong("SAS map.")),
      em("X axis or structural similarity."),
      p("For the SAS maps you can choose the fingerprint you want to use (ECFP, PubChem, MACCS keys).
        These fingerprints are calculated using the R package called rcdk. More information on this package and
        fingerprints can be found", style="text-align:justify", 
        a("here.", href="https://cran.r-project.org/web/packages/rcdk/rcdk.pdf", target="_blank"), 
        "For topological fingerprints, like ECFP, you can choose if you want to use a diameter of 4 or 6. 
        The structural similarity plotted on the X axis is computed using the fingerprint of your choice and Tanimoto similarity.", 
        style="text-align:justify"),
      em("Y axis or activity difference."),
      p("The activity difference plotted on the Y axis 
        is the absolute difference between the activity values of each pair of compounds.",
        style="text-align:justify"),
      em("Thresholds and areas of the plot."),
      p("You can change the X and Y axis thresholds. Depending on tese thresholds the number of
        compounds that you will download from each area in the plot will change.", style="text-align:justify"),
      p("Your SAS map will be divided in 4 areas:"),
      img(src="Figure_2.png", height=300, width=300,
          style="display: block; margin-left: auto; margin-right: auto"),
      p(span("1.", style="color:blue"),"The bottom left area contains pairs of compounds that have
        low molecular similarity and low activity difference, i.e., similarity cliffs.", style="text-align:justify"),
      p(span("2.", style="color:blue"), "In the bottom right area, or smooth SAR, you can find pairs of
        compounds that have high molecular similarity and low activity difference.", style="text-align:justify"),
      p(span("3.", style="color:blue"), "The top right area contains activity cliffs i.e., 
        pairs of compounds with high molecular similarity and high activity difference."),
      p(span("4.", style="color:blue"),"The top left area or non-descriptive, contains pairs of compounds that
have low structural similarity 
        and high activity difference.", style="text-align:justify"),
      p("You can download all the data generated to plot the SAS map and the information for all the pairs of 
        compounds on each area of the plot, as well as their similarity, the most active compound in the pair,
        activity difference and SALI", style="text-align:justify", strong("If you change a certain setting you 
                                                                          will have to wait until the new
                                                                          is produced before you can download the new data.")),
      p("If you select a dot or group of dots on the plot, it will tell you which pair or pairs of
        compounds are in that area.", style="text-align:justify"),
      p(style="text-align:justify", "The color of the plot changes depending on the values calculated for the Structure-Activity Landscape Index (SALI)
        and the most active compound in the pair. You can find more information on SALI", style="text-align:justify",a("here.", href="http://pubs.acs.org/doi/abs/10.1021/ci7004093", target="_blank"), 
        "Pairs of compounds with the highest SALI will be colored green, pairs of compounds with intermediate SALI values
        will be orange to yellow and pairs of compounds in red will have the lowest SALI values. Similarly to SALI the most active 
        compounds on each pair will be red, intermediate compounds will be yellow-to-orange and the less active will be green.
        You can find some examples of SAS maps color-coded by the most active compound in the pair", style="text-align:justify", 
        a("here.", href="https://www.ncbi.nlm.nih.gov/pubmed/19434846", style="text-align:justify", target="_blank")),
       h4(strong("DAD maps.")),
      em("X and Y axis Activity Difference 1 and Activity Difference 2."),
      p("For the DAD map you can choose two biological activities from your input file to be plotted. Each point on the 
        plot corresponds to the pairwise activity difference against target one (activity difference 1) and the pairwise
        activity difference agains target 2 (activity difference 2) for each possible pair in the data set.", style="text-align:justify"),
      em("Thresholds and areas on the plot."),
      p("You can change the X and Y axis superior and inferior thresholds. Depending on these thresholds
        the number of compounds that you will download from each area in the plot will change."),
      em("Your DAD map will be divided in 9 areas:"),
      img(src="Figure_4.png", height=300, width=300,
          style="display: block; margin-left: auto; margin-right: auto"),
      p(span("1.", style="color:blue"), "Areas Z1 up and Z1 down or Z1u and Z1d, respectively, 
contain pairs of compounds for which structural changes have a similar impact on the activity towards the two targets. 
For Z1d the structural changes decrease 
        the activity for the two targets and for Z1u the structural changes increase the activity for the two targets.", style="text-align:justify"),
      p(span("2.", style="color:blue"), "Areas Z2 up and Z2 down or Z2u and Z2d, respectively: 
        indicate that the change in activity for the compounds in the pair is opposite for each target,
        i.e., changes in structure increase the activity for one target, while decreasing activity for the other
        target. For Z2u, structural changes increased the activity on the second target and decreased the activity
        in the first target. Z2d, structural changes decrease the activity for the second target and increase the 
        activity for the first target.", style="text-align:justify"),
      p(span("3.", style="color:blue"), "Areas Z3 and Z4 or Z3d, Z3u and Z4l, Z4r: indicate that structural changes result 
        in significant changes in activity towards one target, but not a large change towards the other target.", style="text-align:justify"),
      p(span("4.", style="color:blue"),"4.	Area Z5 indicates pair of compounds with similar activity against target one and target 2, therefore, 
        structural changes have little or no impact on the activity against the two targets. Compounds in Z5 might be of significant interest if 
        the goal of the analysis is to identify promiscuous hits.", style="text-align:justify"),
      p("The classification of data points in an activity-difference map is independent of the structure similarity. However, in the DAD map the 
        structural similarity of each pair of compounds can be mapped using a continuous color scale: compounds with high similarity will be colored green,
        compounds with intermediate similarity will be colored yellow to orange and compounds with low similarity will be red. You can also add color to the
        data points with the selectivity option. For this option, if a data point is colored red this indicates that one of the compounds in that pair could be
        selective to one of the targets; yellow to orange, one of compounds is moderately selective; green, there are no important activity differences.", style="text-align:justify"),
      p("You can download all the data used to plot the DAD map and the information for all the pairs of compounds on each area of the plot, as well as their selectivity, similarity, 
        activity difference 1 and activity difference 2.", style="text-align:justify", strong("If you change a certain setting you will have to wait until the new plot is produced 
                                                                                              before you can download the new data.", style="text-align:justify")),
      p("If you select a data point or a group of data points on the plot, it will tell you which pair or pairs of compounds are in that area."),
      p("You can find more information about SAS and DAD maps", a("here.", href="http://pubs.acs.org/doi/pdf/10.1021/ci300362x", target="_blank"))
      
    )
      
    ),
    tabPanel("Template",
             sidebarPanel(
               p("Here you can download the .csv file with the
                 columns in the order the App will need, you can fill the template with your data.", style="text-align:justify"),
             p("If your .csv is saved as separated with semicolon (;), it will give an error.
               Make sure it is saved as separated with comma (,).",  style="text-align:justify"),
             p("Click", a("here", href="https://www.difacquim.com/d-tools/", target="_blank"),
               "to download the tamplate."),
             p("You can also download", strong("two example data sets"), "to test this App. The example data sets are published",
               a("here.", href="http://pubs.rsc.org/en/content/articlelanding/2016/ra/c6ra07224k#!divAbstract",target="_blank")),
             p("Click", a("here", href="https://www.difacquim.com/d-tools/", target="_blank"), 
               "to download the examples.", style="text-align:justify")
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
                  choices = c("SALI", "Max.Activity","Density")),
      numericInput("ablineX", "Choose the X axis threshold",
                   min = 0, max = 1, value =0.5, step = 0.5),
      numericInput("ablineY", "Choose the Y axis threshold",
                   min = 0, max = 5, value = 2, step = 0.5),
      downloadButton('downloadData', 'Download SASmap data')
    
    ),
    mainPanel(
      plotOutput("sasplot", width = 900, height = 700,  brush = brushOpts(id = "plot_brush")),
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
                           choices = c("Similarity", "Selectivity","Density")),
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
  tabPanel("Contact",
           sidebarPanel(
                        br(),
                        p("We would appreciate your feedback!"),
                        p("We want to keep improving all out Apps to make them easier to use, 
                          if you have any suggestions or if you experience an error, please send
                          an email to:",style="text-align:justify"),
                        p(a("Mariana Gonzalez Medina", href="https://www.researchgate.net/profile/Mariana_Gonzalez-Medina2",
                            target="_blank"), strong("mgm_14392@comunidad.unam.mx")),
                        p(a("Oscar Mendez Lucio", href="https://www.researchgate.net/profile/Oscar_Mendez-Lucio", 
                            target="_blank"), strong("oscarmen@comunidad.unam.mx")),
                        p(a("Jose L. Medina Franco",href="https://www.researchgate.net/profile/Jose_Medina-Franco/publications",
                            target="_blank"), strong("medinajl@unam.mx"))
                        
           )
  ),
  tabPanel("Acknowledgements",
           sidebarPanel(
                        h3("Funding"),
                        p("UNAM: PAPIME PE200116; PAIP 5000-9163")
                        
           )
  ),
tabPanel("D-Tools",
         h3("DIFACQUIM tools for Cheminformatics"),
         fluidRow(
           column(6,
                  a(img(src="CDPs.png", height = 300, width = 450),
                    href="https://consensusdiversityplots-difacquim-unam.shinyapps.io/RscriptsCDPlots/",
                    target="_blank")),
           column(6,
                  a(img(src="UGFig5.png", height = 300, width = 400),
                    href="http://132.248.103.152:3838/PUMA/",
                    target="_blank"))
         ),
         fluidRow(
           column(6,
                  a("Consensus Diversity Plots (CDPs)",
                    href="https://consensusdiversityplots-difacquim-unam.shinyapps.io/RscriptsCDPlots/",
                    target="_blank")),
           column(6,
                  a("Platform for Unified Molecular Analysis (PUMA)",
                    herf="http://132.248.103.152:3838/PUMA/",
                    target="_blank")))
)
  
  )
)
)
