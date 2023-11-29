
source(paste(getwd(),'global2v.r',sep="/"))


setSliderColor <- function(color, sliderId) {
  
  # some tests to control inputs
  stopifnot(!is.null(color))
  stopifnot(is.character(color))
  stopifnot(is.numeric(sliderId))
  stopifnot(!is.null(sliderId))
  
  
  sliderId <- sliderId - 1
  
  # create custom css background for each slider
  # selected by the user
  sliderCol <- lapply(sliderId, FUN = function(i) {
    paste0(
      ".js-irs-", i, " .irs-single,",
      " .js-irs-", i, " .irs-from,",
      " .js-irs-", i, " .irs-to,",
      " .js-irs-", i, " .irs-bar-edge,",
      " .js-irs-", i,
      " .irs-bar{  border-color: transparent;background: ", color[i+1],
      "; border-top: 1px solid ", color[i+1],
      "; border-bottom: 1px solid ", color[i+1],
      ";}"
    )
  })
  
 
  custom_head <- tags$head(tags$style(HTML(as.character(sliderCol))))
  return(custom_head)
}

sidebar <- dashboardSidebar(
  sidebarMenu(
   
    menuItem("INFO", icon =  tags$i(  class="fa-solid fa-house fa-lg"),tabName = "dashboard"),
    menuItem("Data preprocessing", icon = tags$i(class="fa-solid fa-cloud-arrow-up fa-lg"), tabName = "widgets1"),
    menuItem("Differential expression analysis", icon = tags$i(class="fa-solid fa-magnifying-glass-chart fa-lg"), tabName = "charts", startExpanded = TRUE),
    menuItem("Plots of Differential expression ",icon =   tags$i(  class="fa-solid fa-chart-line fa-lg"), tabName='charts2'),
    menuItem("Gene Set enrichment analysis ", icon = tags$i(  class="fa-solid fa-magnifying-glass-chart fa-lg"), tabName = "widgets2")
  )
)



body <- dashboardBody(
  
  
  
  
  
  tabItems(
    
    
    tabItem(tabName = "dashboard",
            fluidRow(
              
              
              box(title = "User Guide", width = 12, solidHeader = T, status = "info",
                  column(10,
                         includeMarkdown("Shiny_App_Manual.RMD")
                  )
              )
              
              
            )) ,
    tabItem(tabName = "widgets1",
            ############################################# Trimming Box ################################
            box(title = "Raw Data Preprocessing", solidHeader = T, status = "primary", width = 12,
                hr(),
                fluidRow(
                  add_busy_gif(
                    src="https://i.gifer.com/4vsm.gif",
                    timeout = 100,
                    position = c( "full-page"),
                    margins = c(10, 10),
                    overlay_color = "rgba(0, 0, 0, 0.5)",
                    overlay_css = NULL,
                    height = "50px",
                    width = "50px"
                  ),
                  
                  
                  
                  tags$head(
                    tags$style(
                      HTML(".shiny-notification {
                  
                 
                 

                 background-color:#bfb7ee;
             }"
                      )
                    )
                  ),
                  
                  
                  
                  column(6,
                         
                         box(title = "Upload files ", solidHeader = T, status = "primary", width = 10, collapsible = T,id = "inputbox",
                             h5(strong("First you need to upload the files that you want to use in the trimming process")),
                             fileInput('file1', 'Choose  Files',multiple = TRUE,
                                       accept=c('.fq','.fa',
                                                '.fastq' , '.fasta')) ,
                             
                             withSpinner( verbatimTextOutput('result')),
                             h5(strong("Optional:")),
                             actionButton("delete", "Delete all files"),
                         )
                  ),
                  column(6,
                         box(title = "Trimming process", solidHeader = T, status = "primary", width = 10, collapsible = T,id = "trimmingbox",
                             
                             
                             selectizeInput("pair1", label = "Choose R1/Forward reads file",choices = NULL ),
                             
                             selectizeInput("pair2", label = "Choose R2/Reverse reads file",choices = NULL),
                             
                             
                             actionButton("initRbowtie","Remove adaptors", class = "btn-info", style = "width: 100%")
                             
                             
                             
                         )
                         
                  ),
                )
            )
    ),
    
    tabItem(tabName = "charts",
            
            
            
            fluidRow(
              
              add_busy_gif(
                src="https://i.gifer.com/4vsm.gif",
                timeout = 100,
                position = c( "full-page"),
                margins = c(10, 10),
                overlay_color = "rgba(0, 0, 0, 0.5)",
                overlay_css = NULL,
                height = "50px",
                width = "50px"
              ),
              
              
              
              column(6,
                     ############################################### Mapping Box ##################################
                     box(title = "Mapping process", solidHeader = T, status = "primary", width = 12, collapsible = T,id = "mappingbox",
                         
                         
                         
                         box(title = "Build index", solidHeader = T, status = "primary", width = 7, collapsible = T,id = "indexbox",
                             h5(strong("Build the genome index if you haven't already done it.")),
                             
                             
                             selectizeInput("reference_fasta", label = "Choose reference genome (.fasta/.fna/.fa):",choices = NULL ),
                             
                             actionButton("initBuild","Build index", class = "btn-info", style = "width: 100%")),
                         
                         
                         box(title = "Alignment", solidHeader = T, status = "primary", width = 10, collapsible = T,id = "alignmentbox",     
                             wellPanel(
                               
                               
                               
                               column(10,
                                      numericInput("minFrag","Select the minimum length of the fragment (minFragLength):", value = 50)
                               ),
                               
                               column(10,
                                      numericInput("maxFrag","Select the maximum length of the fragment (maxFragLength):", value = 600)
                               ),
                               column(10,
                                      numericInput("misMatch","Select the maximum number of mis-matched bases allowed in the alignment (maxMismatches):", value = 3)
                               ),
                               
                               column(10,
                                      selectInput("format","Select the output format (OutputFormat):", choices = c("SAM", "BAM"), selected = "BAM")
                               ),
                               column(10,
                                      numericInput("nthreads","Select the number of the threads you want to use, depending on your resources (nthreads):", value = 12, max = 40)
                               ),
                               
                               #h5(strong("Proceed to alignment")),
                               
                               actionButton("initAlign","Alignment", class = "btn-info", style = "width: 100%")
                               
                               
                               
                             )
                             )
                     )),
              
              ############################################### Quantification Box ##################################
              
              column(6,
                     box(title = "Quantification process", solidHeader = T, status = "primary", width = 12, collapsible = T,id = "quantyfbox",
                         
                         selectizeInput("reference_gtf", label = "Choose reference genome (.gtf):",choices = NULL ),
                         
                         wellPanel(
                           #column(4,
                           #selectInput("isGTFAnnotationFile","isGTFAnnotationFile:", choices = c(TRUE, FALSE), selected = TRUE)
                           # ),
                           
                           column(10,
                                  selectInput("overlap","Is a read allowed to be assigned to more than one feature, if an overlap is found? (AllowMultiOverlap):", choices = c(TRUE, FALSE), selected = FALSE)
                           ),
                           
                           
                           column(10,
                                  numericInput("minOverlap","Select the minimum number of overlapped bases required for assigning a read to a feature (minOverlap):", value = 1)
                           ),
                           column(10,
                                  selectInput("countMultiMappingReads","Should multi-mapping reads/fragments be counted? (countMultiMappingReads):", choices = c(TRUE, FALSE), selected = TRUE)
                           ),
                           
                           
                           
                           column(10,
                                  numericInput("minFragL","Select the minimum length of the fragment (minFragLength):", value = 50)
                           ),
                           
                           column(10,
                                  numericInput("maxFragL","Select the maximum length of the fragment (maxFragLength):", value = 600)
                           ),
                           column(10,
                                  numericInput("nthreads","Select the number of the threads you want to use, depending on your resources (nthreads):",value = 12, max = 40)
                           ),
                           
                           
                           
                           
                           column(10,
                           actionButton("initQuantify","Quantification", class = "btn-info", style = "width: 100%")),
                           
                           downloadButton('downloadcountCSV','Save Results as CSV File', class = "btn btn-info", style="margin: 7px;")
                           
                           
                         ))
                     
              ),  
              
              column(10,
                     tags$div(class = "BoxArea2",
                              br(),
                              withSpinner(dataTableOutput('Table')),
                              br(),
                     ),
                     tags$div(class = "clearBoth")
              ),
              tags$div(class = "clearBoth"),
              
              
             
              
             
             ############################################### Differential Expression Analysis Box ##################################          
           
               box(title = "Differential Expression Analysis", solidHeader = T, status = "success", width = 12,
                  hr(),   
                  column(8,
                         box(title = "Upload the counts file", solidHeader = T, status = "success", width = 8, collapsible = T,id = "uploadcsvbox",
                             h5(strong("Make sure that the header of the csv contains only the  sample names")),                
                             
                             
                             fileInput('target_upload_csv', 'Choose file to upload',
                                       accept = c(
                                         'text/csv',
                                         'text/comma-separated-values',
                                         '.csv'
                                       )),
                             
                             
                             
                         ),
                         br(),
                         column(8,
                                br(),
                         DT::dataTableOutput("count_table"))),
                         
                
                  column(8,
                         br(),
                         box(title = "OPTION 1: Upload the phenodata file ", solidHeader = T, status = "success", width = 8, collapsible = T,id = "uploadcsvbox2",
                             h5(strong("Make sure that the file contains the same samples as the counts file and that the names match")),
                             
                             fileInput('metadatafile', '',
                                       accept=c('text/csv', 
                                                'text/comma-separated-values,text/plain', 
                                                '.csv'),multiple = FALSE
                             )
                         )),
                  
                  
                 
                                        
                                           
                                       column(8,
                                              box(title = "OPTION 2: Create phenodata table", solidHeader = T, status = "success", width = 10, collapsible = T,id = "phenodatatable",
                                                  setSliderColor("#28c5cd", 1),    
                                                  column(8,
                                                         
                                                  radioButtons("n", "Choose the number of columns for the dataframe:", c("1", "2",'3'))
                                                  #sliderInput("n",label= "Choose the number of columns for the dataframe:", min=1, max=3, value=1, ticks=TRUE)
                                                  ,
                                             
                                                       h5(strong("Make sure that the samples and the conditions are in the right order.")),
                                                       #h4("Add Conditions/Factors"),
                                                       textInput("conditionName1", "Column Name 1:", placeholder = "Eg. Condition"),
                                                       textInput(  "conditions1",  "List of Conditions/Factors (comma seperated):", placeholder = "Eg. control, treatment" ),
                                                       textInput("conditionName2", "Column Name 2:", placeholder = "Eg. Time:"),
                                                       textInput(  "conditions2",  "List of Conditions/Factors (comma seperated):", placeholder = "Eg. 1hr,2hr" ),
                                                       textInput("conditionName3", "Column Name 3:", placeholder = "Eg. Sex"),
                                                       textInput(  "conditions3",  "List of Conditions/Factors (comma seperated):", placeholder = "Eg. male, female" ),
                                                        br(), 
                                                        actionButton(
                                                         "addConditions",
                                                         "Create table",
                                                         class = "btn btn-primary"),
                                                        br(),
                                                ),
                                                
                                               
                                     br(),
                                     column(8,
                                            br(),
                                     DT::dataTableOutput("phenodata_table")))),
                                     
              
                  
                                          
                                            
                   
                  
                 
                  
                  
                  column(8,
                         box(title = "DESeq2 Object Settings  ", solidHeader = T, status = "success", width = 10, collapsible = T,id = "DDSBox",
                            
                             wellPanel(
                             
                             h5(strong("Design your own formula depending on your analysis needs (differential expression with multiple or single factors)")), 
                             
                             textInput("designFormula","Design Formula:", placeholder = "~ Conditions"),
                             
                             
                            
                             h5(strong(" You can specify the comparison of interest, and the levels to compare.
                                           The level given last is the base level for the comparison.")), 
                             
                             selectInput("factorNameInput",label ="Choose the contrast column:", choices = NULL),
                             actionButton("choice", "Update Levels by your contrast choice"),
                             
                             selectInput("condition1",label ="Choose Level 1:", choices = NULL),
                             
                             selectInput("condition2",label ="Choose Level 2:", choices = NULL),
                             column(4,
                                    numericInput("padj","P-Adjust Cutoff:", value = 0.05)
                             ),
                             column(4,
                                    numericInput("log2F","log2FoldChange Cutoff:", value = 1)
                             ),
                             actionButton("initResults","Get Results", class = "btn-info", style = "width: 100%"),
                             downloadButton('downloadDDSCSV2','Save Results as CSV File', class = "btn btn-info", style="margin: 7px;")),
  
                               
                            
                           
                         )),
                  
                  column(10,
                        
                       DT::dataTableOutput("Table2")),
                  br(),
                  
                        
                  )
                  
              )
              
            ),
    
    
    ####################################### Differential Expression Plots Box #####################################
    tabItem(tabName = "charts2",  
            
            fluidRow( 
              
              add_busy_gif(
                src="https://i.gifer.com/4vsm.gif",
                timeout = 100,
                position = c( "full-page"),
                margins = c(10, 10),
                overlay_color = "rgba(0, 0, 0, 0.5)",
                overlay_css = NULL,
                height = "50px",
                width = "50px"
              ),
              
              
              box(title = "Color palette", solidHeader = T, status = "success", width = 8, collapsible = T,id = "palette",
                  
                  column(4,
                         colourInput("col1", "Select color 1", "#30ABAD", allowTransparent = T)),
                  column(4,
                         colourInput("col2", "Select color 2", "#BF58DB", allowTransparent = T)),
                  column(4,
                         colourInput("col3", "Select color 3", "#555C5C", allowTransparent = T))),
              
              
              box(title = "Volcano plot", solidHeader = T, status = "success", width = 10, collapsible = T,id = "volcanoplot2",
                  
                  column(9,
                         wellPanel(
                           withSpinner(plotOutput(outputId = "volcanoplot2"))
                         )
                  ),
                  downloadLink("downloadPlotvol", "Download Plot") 
              ),
              
              box(title = "Ma plot", solidHeader = T, status = "success", width = 10, collapsible = T,id = "maplot2",
                  
                  column(9,
                         wellPanel(
                           withSpinner(plotOutput(outputId = "maplot2"))
                         )
                  ),
                  downloadLink("downloadPlotma", "Download Plot") 
              ),
              
              box(title = "PCA plot", solidHeader = T, status = "success", width = 10, collapsible = T,id = "pcaplot",
                  fluidRow( 
                    column(6,
                           h4(strong("Plot Settings:")),
                           wellPanel(
                             
                             fluidRow(
                               column(12,
                                      selectInput("group","Group of interest:", choices = c())
                               )))),
                    
                    box(title = "PCA ", solidHeader = T, status = "success", width = 10, collapsible = T,id = "pcaplot2",          
                        column(10,
                               wellPanel(
                                 withSpinner(plotOutput(outputId = "pcaplot"))
                               ),
                               
                        ),
                        
                        downloadLink("downloadPlot7", "Download Plot"))
                  )),
              
              box(title = "Heatmap ", solidHeader = T, status = "success", width = 10, collapsible = T,id = "heatmap",
                  fluidRow( 
                    column(6,
                           h4(strong("Plot Settings:")),
                           wellPanel(
                             
                             fluidRow(
                               
                               column(12,
                                      selectInput("group1","Group 1:", choices = c())
                               ),
                               column(12,
                                      selectInput("group2","Group 2 (Optional):", choices = c())
                               ),
                               column(12,
                                      radioButtons("m", "Number of groups:", c("1", "2")))
                             ),
                             
                             
                             
                             div(style = "clear:both;")
                           )
                    ),
                    
                    
                    box(title = "Heatmap", solidHeader = T, status = "success", width = 11, collapsible = T,id = "heatmap",
                        
                        column(10,
                               wellPanel(
                                 withSpinner(plotOutput(outputId = "heatmap"))
                               )
                        ),
                        downloadLink("downloadPlot6", "Download Plot") 
                    )
                  ),
                  
                  
              ),
              box(title = "Gene Expression Boxplot", solidHeader = T, status = "success", width = 12,
                  
                  fluidRow(
                    column(6,
                           wellPanel(
                             column(12, 
                                    selectizeInput("sel_gene",
                                                   label="Gene Name/Id (Select 1 or more)",
                                                   choices = NULL,
                                                   multiple=TRUE,
                                                   options = list(
                                                     placeholder = 
                                                       'Type to search for a gene name/id'
                                                   ) #,
                                    )
                             ),
                             
                             
                             
                             div(style = "clear:both;")
                           )
                    ),
                    column(6,
                           h4(strong("Plot Settings:")),
                           wellPanel(
                             
                             fluidRow(
                               column(12,
                                      selectInput("boxplotX","X-axis", choices = c())
                               ),
                               column(12,
                                      selectInput("boxplotFill", "Grouped by", choices = c())
                               )
                             ),
                             
                             div(style = "clear:both;")
                           )
                    ),
                    
                    
                    
                    box(title = "Boxplot", solidHeader = T, status = "success", width = 11, collapsible = T,id = "boxPlot",
                        
                        column(10,
                               wellPanel(
                                 withSpinner(plotOutput(outputId = "boxPlot"))
                               )
                        ),
                        downloadLink("downloadPlot8", "Download Plot")
                    ))),
              
              box(title = "Count Boxplot", solidHeader = T, status = "success", width = 11, collapsible = T,id = "boxplot2",
                  column(6,
                         h4(strong("Choose color palette:")),
                         wellPanel(
                           
                           fluidRow(
                             column(12,
                                    selectInput("Palette","Palette:", c("BuGn", "Blues", "OrRd","Spectral","Purples","YlGnBu","YlOrBr"), selected = "BuGn")
                             )))),
                  column(10,
                         wellPanel(
                           withSpinner(plotOutput(outputId = "boxplot2"))
                         )
                  ),
                  
              )
              
              
              
              
            )
            
    ),
    
    
    
    tabItem(tabName = "widgets2",
            
            
            fluidRow(
              
              add_busy_gif(
                src="https://i.gifer.com/4vsm.gif",
                timeout = 100,
                position = c( "full-page"),
                margins = c(10, 10),
                overlay_color = "rgba(0, 0, 0, 0.5)",
                overlay_css = NULL,
                height = "50px",
                width = "50px"
              ),
              
              #box(title = "Gene Enrichement Analysis", solidHeader = T, status = "success", width = 12,
              # hr(),
              #column(8,
                     box(title = "Upload Data", solidHeader = T, status = "success", width = 8, collapsible = T,id = "uploadbox",
                         
                         h5(strong("Before you upload your dataset make sure that the column that contains the entrezids
                   is named  gene_name and that the column that contains the log2 Fold Change results is named log2FoldChange.")), 
                         
                         
                         
                         fileInput('target_upload3', 'Choose file to upload',
                                   accept = c(
                                     'text/csv',
                                     'text/comma-separated-values',
                                     '.csv'
                                   )),
                         
                         h5(strong("In case there is an  pattern in front of the gene ids in the gene_name column, you should remove it.")), 
                         
                         
                         column(4,
                                textInput("pattern",label="Give the  pattern: ")),
                         column(10,
                                checkboxInput("removePattern","Remove pattern", value = F))
                     ),
                     
                     column(8,
                     DT::dataTableOutput("sample_table"),
                     br(),
                     ),
              
              ############################################# Annotation Box ##############################            
              
              box(title = "Annotation  Parameters", solidHeader = T, status = "primary", width = 8, collapsible = T,id = "Annotationhub",
                  
                  fluidRow(
                    column(10,
                           textInput("organism",label="Select organism (latin name):")),
                    
                    column(10,
                           actionButton("search","Search in Annotationhub", class = "btn-info")),
                    column(10,
                           withSpinner( verbatimTextOutput('text'))),
                    column(10,
                           textInput("annot",label="Select annotation file code:")),
                    column(10,
                           actionButton("annotation","Download file from Annotationhub", class = "btn-info")),
                  )),
              
              
              
              
              
              
              ###################################### Enrichment Box ###################################              
              box(title = "Gene Enrichement Object Parameters", solidHeader = T, status = "primary", width = 8, collapsible = T,id = "createGoBox",
                  
                  wellPanel(
                    # column(4,
                    #  selectizeInput("keytypes","Keytypes:", choices = c("MF", "BP", "CC","ALL"), selected = "ALL")
                    #  ),
                    column(10,
                           
                           selectizeInput("org",
                                          label="Select organism (short name):",
                                          choices = NULL,
                                          
                                          options = list(
                                            placeholder = 
                                              'Type to search for an organism 
                                     ex:human'
                                          )) ),
                    
                    
                    column(10,
                           selectInput("ontology","Select Ontology:", choices = c("MF", "BP", "CC","ALL"), selected = "ALL")
                    ),
                   # column(10,
                           #numericInput("nPerm","Test the significance of gene set enrichment, by increasing the permutation parameter (Permutation) :", value = 1000, min = 1, max = 100000)
                  #  ),
                    column(10,
                           numericInput("minGSSize","Select the minimal size of each geneSet for analyzing (minGSSize):", value = 3)
                    ),
                    column(10,
                           numericInput("maxGSSize","Select maximal size of genes annotated for testing (maxGSSize):", value = 800)
                    ),
                    column(10,
                           numericInput("pvalCuttoff","P-Value Cutoff:", value = 0.05)
                    ),
                    column(10,
                           selectInput("pAdjustMethod","Select the p-value adjustment method (pAdjustMethod):", choices = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"), selected = "none")
                    ),
                    
                   
                  column(10,
                  actionButton("initGo","Create gseGO Object", class = "btn-info", style = "width: 100%")),
                 
                  downloadButton('downloadgseGoCSV','Save Results as CSV File', class = "btn btn-info", style="margin: 7px;")
                  ),
                  
                  
                  
                  fluidRow(
                    
                    
                    column(10,
                           tags$div(class = "BoxArea2",
                                    withSpinner(dataTableOutput('gseGoTable'))
                           ),
                           tags$div(class = "clearBoth")
                    ),
                    tags$div(class = "clearBoth")
                    
                  ), 
              ),
              
              
              
              
              ###################################### Enrichment Plots Box ####################################    
              
              box(title = "Dot Plot", solidHeader = T, status = "success", width = 10, collapsible = T,id = "dotPlot",
                  fluidRow(
                    column(3,
                           wellPanel(
                             numericInput("showCategory_dot", "number of categories to show", value = 10)
                           )
                    ),
                    column(6,
                           wellPanel(
                             withSpinner(plotOutput(outputId = "dotPlot"),type = 8)
                             
                           ),
                           
                           
                           
                           
                    ),
                    downloadLink("downloadPlot", "Download Plot")
                  )
              ),
              
              box(title = "Category Netplot", solidHeader = T, status = "success", width = 10, collapsible = T,id = "cnetplot",
                  fluidRow(
                    column(9,
                           wellPanel(
                             numericInput("showCategory_cne", "number of categories to show", value = 5)
                           )
                    )
                    ,
                    column(9,
                           wellPanel(
                             plotOutput(outputId = "cnetplot")
                           )
                    ),
                    downloadLink("downloadPlot4", "Download Plot")
                  )
              ),
              
              box(title = "Heat Plot", solidHeader = T, status = "success", width = 10, collapsible = T,id = "heatplot",
                  
                  column(9,
                         wellPanel(
                           withSpinner(plotOutput(outputId = "heatplot"))
                         )
                  ),
                  downloadLink("downloadPlot5", "Download Plot")
                  
              ),
              
              
              box(title = "Ridge Plot", solidHeader = T, status = "success", width = 10, collapsible = T,id = "ridgeplot",
                  fluidRow(
                    column(3,
                           wellPanel(
                             numericInput("showCategory_ridge", "number of categories to show", value = 10)
                           )
                    ),
                    column(6,
                           wellPanel(
                             withSpinner(plotOutput(outputId = "ridgePlot"),type = 8)
                           )
                    ),
                    downloadLink("downloadPlot2", "Download Plot")
                  )
              ),
              box(title = "GSEA Plot", solidHeader = T, status = "success", width = 10, collapsible = T,id = "gseaplot",
                  fluidRow(
                    column(9,
                           wellPanel(
                             numericInput("geneSetId_gsea", "Gene Set ID", value = 1)
                           )
                    ),
                    column(9,
                           wellPanel(
                             plotOutput(outputId = "gseaplot",width = "100%", height = "400px")
                           )
                    ),
                    downloadLink("downloadPlot3", "Download Plot")
                  )
              )
            )
    )
    
    
    
  ),
  
  
  customTheme <- shinyDashboardThemeDIY(
    
    
    appFontFamily = "Arial"
    ,appFontColor = "rgb(0,0,0)"
    ,primaryFontColor = "rgb(0,0,0)"
    ,infoFontColor = "rgb(0,0,0)"
    ,successFontColor = "rgb(0,0,0)"
    ,warningFontColor = "rgb(0,0,0)"
    ,dangerFontColor = "rgb(0,0,0)"
    ,bodyBackColor = "rgb(244, 246, 246 )"
    
    # header
    ,logoBackColor = "rgb(255,255,255)"
    
    ,headerButtonBackColor = "rgb(255,255,255 )"
    ,headerButtonIconColor = "rgb(40,	197	,205)"
    ,headerButtonBackColorHover = "rgb(127, 255, 212)"
    ,headerButtonIconColorHover = "rgb(224,224,224)"
    
    ,headerBackColor = "rgb(255,255,255 )"
    ,headerBoxShadowColor = "#aaaaaa"
    ,headerBoxShadowSize = "2px 2px 2px"
    
    # sidebar
    ,sidebarBackColor = "rgb(255,255,255 )"
    ,sidebarPadding = 0
    
    ,sidebarMenuBackColor = "rgb(255,255,255 )"
    ,sidebarMenuPadding = 0
    ,sidebarMenuBorderRadius = 0
    
    ,sidebarShadowRadius = "3px 8px 9px"
    ,sidebarShadowColor = "#aaaaaa"
    
    ,sidebarUserTextColor = "rgb(40,	197	,205)"
    ,sidebarSearchBackColor = "rgb(23,103,124)"
    ,sidebarSearchIconColor = "rgb(255,255,255)"
    ,sidebarSearchBorderColor = "rgb(255,255,255)"
    
    ,sidebarTabTextColor = "rgb(40,	197,205)"
    ,sidebarTabTextSize = 13
    ,sidebarTabBorderStyle = "none none solid none"
    ,sidebarTabBorderColor = "rgb(224,224,224)"
    ,sidebarTabBorderWidth = 1
    
    ,sidebarTabBackColorSelected = cssGradientThreeColors(
      direction = "right"
      ,colorStart = "rgb(127, 255, 212)"
      ,colorMiddle = "rgb(127, 255, 212)"
      ,colorEnd = "rgb(127, 255, 212)"
      ,colorStartPos = 0
      ,colorMiddlePos = 30
      ,colorEndPos = 100
    )
    ,sidebarTabTextColorSelected = "rgb(0,0,0)"
    ,sidebarTabRadiusSelected = "0px 20px 20px 0px"
    
    ,sidebarTabBackColorHover = cssGradientThreeColors(
      direction = "right"
      ,colorStart = "rgb(127, 255, 212)"
      ,colorMiddle = "rgb(127, 255, 212)"
      ,colorEnd = "rgb(127, 255, 212)"
      ,colorStartPos = 0
      ,colorMiddlePos = 30
      ,colorEndPos = 100
    )
    ,sidebarTabTextColorHover = "rgb(50,50,50)"
    ,sidebarTabBorderStyleHover = "none none solid none"
    ,sidebarTabBorderColorHover = "rgb(75,126,151)"
    ,sidebarTabBorderWidthHover = 1
    ,sidebarTabRadiusHover = "0px 20px 20px 0px"
    
    # boxes
    ,boxBackColor = "rgb(255,255,255)"
    ,boxBorderRadius = 12
    ,boxShadowSize = "3px 5px 5px"
    ,boxShadowColor = "#aaaaaa"
    ,boxTitleSize = 16
    ,boxDefaultColor = "rgb(210,214,220)"
    ,boxPrimaryColor = "rgb(40,	197	,205)"
    ,boxInfoColor = "rgb(40,	197	,205)"
    ,boxSuccessColor = "rgb(40,	197	,205)"
    ,boxWarningColor = "rgb(244,156,104)"
    ,boxDangerColor = "rgb(255,88,55)"
    
    ,tabBoxTabColor = "rgb(255,255,255)"
    ,tabBoxTabTextSize = 14
    ,tabBoxTabTextColor = "rgb(0,0,0)"
    ,tabBoxTabTextColorSelected = "rgb(0,0,0)"
    ,tabBoxBackColor = "rgb(255,255,255)"
    ,tabBoxHighlightColor = "rgb(127, 255, 212)"
    ,tabBoxBorderRadius = 0
    
    # inputs
    ,buttonBackColor = "rgb(127, 255, 212)"
    ,buttonTextColor = "rgb(0,0,0)"
    ,buttonBorderColor = "rgb(200,200,200)"
    ,buttonBorderRadius = 15
    
    ,buttonBackColorHover = "rgb(56,161,187)"
    ,buttonTextColorHover = "rgb(127, 255, 212)"
    ,buttonBorderColorHover = "rgb(56,161,187)"
    
    ,textboxBackColor = "rgb(255,255,255)"
    ,textboxBorderColor = "rgb(200,200,200)"
    ,textboxBorderRadius = 15
    ,textboxBackColorSelect = "rgb(245,245,245)"
    ,textboxBorderColorSelect = "rgb(200,200,200)"
    
    # tables
    ,tableBackColor = "rgb(255,255,255)"
    ,tableBorderColor = "rgb(240,240,240)"
    ,tableBorderTopSize = 1
    ,tableBorderRowSize = 1
    
    
    
  )
  
)


ui <-dashboardPage(
  dashboardHeader(title =tags$a(tags$img(src='simplarLOGO.png',height='50px',width='150px'))),
  sidebar,

  body,
  tags$script(src = "https://kit.fontawesome.com/816b4e6b6f.js")
)


#max file size is 12 GB (for the uploaded files)
options(shiny.maxRequestSize = 12000*1024^2) 

server <- function(input, output, session) {
  
  
  
  
  #path to Rbowtie2/extdata/adrm
  destDir <- system.file("extdata", "adrm",package = "Rbowtie2") 
  
  
  #copy user's raw .fastq files inside the subfolder adrm
  output$result <- renderPrint({
    
    inFile <- input$file1
    if (is.null(inFile)) {
      cat("No additional document has been submitted.\n")
      #Update the files inside  adrm and render the names in the textoutput
      result <- list.files(system.file("extdata", "adrm", package="Rbowtie2"),
                           pattern="*.fastq$")
      updateForm()
      return(result)
    }
    cat("Reading file:", inFile$name, "\n")
    cat("size:", inFile$size, " Bytes, type:", inFile$type, "\n")
    if (dir.exists(destDir)){
      cat("Copying file to:", destDir,"\n")
      result <- file.copy( inFile$datapath,
                           file.path(destDir, inFile$name) )
      result<- list.files(system.file("extdata", "adrm", package="Rbowtie2"),
                          pattern="*.fastq$")
      
      
    } else {
      result <- FALSE
    }
    
    updateForm()
    result
  })
  
  observe({
    
    delete_files()
    
  })
  
  #delete all the files from extdata\adrm  after you are done with the trimming process
  delete_files <- eventReactive(input$delete,{
    
    files_to_delete <- dir(path=destDir ,pattern=".fastq")
    file.remove(file.path(destDir, files_to_delete))
   
    showNotification("Files are deleted!", duration = 5,
                     type = ( "message"))
  })
  
  
  
  
  #update the file lists (forward and reverse) in the trimming box in ui
  updateForm = function()
  {
    isolate({
      
      groupv = list.files(system.file(package = "Rbowtie2","extdata","adrm"),pattern = "R1.*.fastq$"  )
      groupv2 = list.files(system.file(package = "Rbowtie2","extdata","adrm"),pattern = "R2.*.fastq$"  )
      
      updateSelectizeInput(session, "pair1", choices = groupv, selected = groupv[1],server=TRUE)
      updateSelectizeInput(session, "pair2", choices = groupv2, selected = groupv2[1],server=TRUE)
      
    }) 
  }
  
  
  
  
  
  ############################# Trimming #################################  
  observe({
    
    Rbowtie2Reactive()
    
  })
  
  td<-getwd()
  
  Rbowtie2Reactive <- eventReactive(input$initRbowtie,{
    withProgress(message = "Trimming , please wait",{ 
      
      pair1<-input$pair1
      pair2<-input$pair2
      reads_1 <- system.file(package = "Rbowtie2","extdata","adrm" , pair1)
      reads_2 <-system.file(package = "Rbowtie2","extdata","adrm" , pair2)
      
      
      (cmdout<-remove_adapters(file1=reads_1,
                               file2=reads_2,
                               adapter1 = "TACACTCTTTCCCTACACGACGCTCTTCCGATCT", 
                               adapter2 ="GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT" ,
                               basename=file.path(td, "trimmed"),
                               
                               overwrite=TRUE,"--threads 3"))
      
      
      #path to "trimmed.pair1.truncated" file
      file1<-file.path(td, "trimmed.pair1.truncated")
      #rename pair1 from "trimmed.pair1.truncated" to "original_file_name_trimmed.fastq" 
      file.rename(file1,(sub("*.fastq$", "_trimmed.fastq", pair1)))
      #path to "trimmed.pair2.truncated" file
      file2<-file.path(td, "trimmed.pair2.truncated")
      #rename pair2 from "trimmed.pair2.truncated" to "original_file_name_trimmed.fastq"
      file.rename(file2,(sub("*.fastq$", "_trimmed.fastq", pair2)))
      #path to "trimmed.settings" file
      file3<-file.path(td, "trimmed.settings")
      file.rename(file3,(sub("*.fastq$", "_trimmed.settings", pair1)))
      
      
      
    })
    showNotification("Trimming completed!", duration = NULL,
                     type = ( "message"))
  })
  
  
  
  
  
  
  
  ############################## Build index ##########################  
  observe({
    
    BuildindexReactive()
    
  })
  #list reference genome file  (.fastq/.fa/.fna)  file from the working directory
  REF_GENOME <- list.files(td, pattern = "*.fasta$|*.fa$|*.fna$" )
  updateSelectizeInput(session, "reference_fasta", choices = REF_GENOME, selected = REF_GENOME[1],server=TRUE)
  BuildindexReactive <- eventReactive(input$initBuild,{
    withProgress(message = "Building index , please wait",{
      
      
      #define the basename for the index
      RSUBREAD_INDEX_BASE <- "genome"  
      
      setProgress(value = 0.3, detail = "This might take a while, please wait ...")
      #build gapped and splitted index to save storage
      buildindex(basename=file.path(td,RSUBREAD_INDEX_BASE),
                 reference=input$reference_fasta,
                 gappedIndex=TRUE,
                 indexSplit=TRUE,
                 memory=4000)   
      
      
      
    })
    showNotification("Completed!", duration = NULL,
                     type = ( "message"))
  })       
  
  
  
  
  
  
  
  ############################# Mapping #################################  
  observe({
    
    MappingReactive()
    
  })
  
  
  MappingReactive <- eventReactive(input$initAlign,{
    
    withProgress(message = "Mapping , please wait",{  
      
      
      # define the basename for the index
      RSUBREAD_INDEX_BASE <- "genome"
      
      #list trimmed forward reads 
      inputfilefwd <-  list.files(td, pattern = "_R1_.*.fastq$" )
      #list trimmed reverse reads 
      inputfilervs <-  list.files(td, pattern = "_R2_.*.fastq$" )
      
      setProgress(value = 0.3, detail = "This might take a while, please wait ...")
      
      align(index=file.path(td,RSUBREAD_INDEX_BASE),
            minFragLength = input$minFrag,
            maxFragLength = input$maxFrag,
            maxMismatches=input$misMatch, 
            nthreads = input$nthreads,
            readfile1=inputfilefwd, 
            readfile2=inputfilervs, 
            output_format=input$format)
      
      
    })
    
    showNotification("Mapping completed!", duration = NULL,
                     type = ( "message"))
    
    
  })
  
  
  
  
  ############################# Quantification #################################
  
  
  #list reference genome files (.gtf) from working directory 
  gtf <- list.files(td, pattern = "*.gtf$" )
  
  updateSelectizeInput(session, "reference_gtf", choices = gtf, selected = gtf[1],server=TRUE)
  
  
  
  
  QuantifyReactive <- eventReactive(input$initQuantify,{
    
    withProgress(message = "Counting , please wait",{ 
      
      
      #list the .SAM or .BAM files from the working directory
      inputfile <- list.files(td, pattern = "*.SAM$|*.BAM$" )
      
      mycounts<-featureCounts(inputfile, 
                              annot.ext= input$reference_gtf,
                              isGTFAnnotationFile=TRUE,
                              isPairedEnd=TRUE,
                              allowMultiOverlap =input$overlap,
                              minOverlap = input$minOverlap,
                              countMultiMappingReads = input$countMultiMappingReads,
                              minFragLength = input$minFragL,
                              maxFragLength = input$maxFragL,
                              nthreads=input$nthreads)
      
      mycounts$counts
      return(list('counts'=mycounts$counts))
      
    })
    
  })
  
  
  
  
  
  #read the counts file
  df_products <- reactive({
    inFile4 <- input$target_upload_csv
    if (is.null(inFile4))
      return(NULL)
    dataframe <- read.csv(inFile4$datapath, header = TRUE,sep = ',')
    
    return(dataframe)
  })
  
  output$count_table<- DT::renderDataTable({
    d <- df_products()
    DT::datatable(d)
  })
  
  

  
  
  
  
  observe({
    tableEditReactive()
  })
  
  tableEditReactive <- eventReactive(input$addConditions,{
    
    if(!is.null(df_products()))
    {
      
      df2<-df_products()
      #for 1 condition
      if (input$n == 1) 
      {
        
        names = colnames(df2)[-1]
        
        column1= input$conditions1
        column1 = strsplit(column1,",")
        
        phenodata <- data.frame(row.names=names,c(column1) )
        colnames(phenodata)[1]  <- input$conditionName1 
        myValues$DF = phenodata

        updateDesignFormula()
        updateFormula()
        return(phenodata)
       }
        
      
        #for 1 condition
        if (input$n == 2) 
        {
          
          
          names = colnames(df2)[-1]
          
          column1= input$conditions1
          column1 = strsplit(column1,",")
          column2= input$conditions2
          column2 = strsplit(column2,",")
          phenodata <- data.frame(row.names=names,c(column1),c(column2))
         
          colnames(phenodata) <- c(input$conditionName1, input$conditionName2)
          myValues$DF = phenodata
          updateDesignFormula()
          updateFormula()
          return(phenodata)
          
          }
 
      
      

    if (input$n == 3) 
    {
      names = colnames(df2)[-1]
      
      column1= input$conditions1
      column1 = strsplit(column1,",")
      column2= input$conditions2
      column2 = strsplit(column2,",")
      column3= input$conditions3
      column3 = strsplit(column3,",")
      phenodata <- data.frame(row.names=names, c(column1),c(column2), c(column3))
      
      colnames(phenodata) <- c(input$conditionName1, input$conditionName2,input$conditionName3)
     
     
      myValues$DF = phenodata
      updateDesignFormula()
      updateFormula()
      return(phenodata)}
      
  
    }
  }
  )
  
  output$phenodata_table<- DT::renderDataTable({
    phenodata=tableEditReactive()
    DT::datatable(phenodata)
  })
  
  
 # iv <- InputValidator$new()
# iv$add_rule("target_upload_csv", sv_required())
 # iv$add_rule("conditionName1", sv_required())
  #iv$add_rule("conditionName1", sv_email())
  #iv$enable()
  
  #read the phenodata file
  observe({
    metadataFileReactive()
  })
  metadataFileReactive <- reactive({
    
    #check if a file is selected, if not then ask to upload a file
    shiny:: validate(
      need( (!is.null(input$metadatafile)),
            message = "Please select a file")
    )
    
    inFile <- input$metadatafile
    if (is.null(inFile))
      return(NULL)
    
    inFile = inFile$datapath  
    
    sep = '\t'
    if(length(inFile) > 0 ){
      testSep = read.csv(inFile[1], header = TRUE, sep = '\t')
      if(ncol(testSep) < 2)
        sep = ','
    }
    else
      return(NULL)
    
    fileContent = read.csv(inFile[1], header = TRUE, sep = sep)
    
    sampleN = colnames(fileContent)[-1]
    metaData <- fileContent[,sampleN]
    metaData <- data.frame(sapply( metaData, as.factor ))
    
    row.names(metaData) <- fileContent[,1]
    if(length(colnames) == 1)
      colnames(metaData)[1] <- sampleN
    myValues$DF = metaData
    
    updateDesignFormula()
    updateFormula()
    
    return(metaData)
  })
  
  
  #Update the design formula in DESeq2 depending on the phenodata file
  updateDesignFormula = function()
  {
    isolate({
      groupvars = colnames(myValues$DF)
      #for 1 condition 
      if(length(groupvars) == 1)
      {
        
        designFormula = paste("~ ",groupvars)
        
      }
      #for more conditions
      else
        designFormula = paste(" ~ ",paste(groupvars, collapse=" + "))
      
      updateTextInput(session,"designFormula", value = designFormula)
      
      
    })
    
  }
  
  #Update contrast choice  
  updateFormula = function()
  {
    isolate({
      
      groupvars = colnames(myValues$DF)
      
      
      updateSelectizeInput(session, "factorNameInput", choices = groupvars, selected = groupvars[1])
      
      
    }) 
  }
  
  
  observe({
    
    Select_name()
    
  })
  
  #Update levels depending on the contrast choice 
  Select_name <- eventReactive(input$choice,{
    
    isolate({
      
      groupvars2 = myValues$DF[,(input$factorNameInput)]
      
      
      updateSelectizeInput(session, "condition1", choices = groupvars2, selected = groupvars2[1])
      updateSelectizeInput(session, "condition2", choices = groupvars2, selected = groupvars2[2])
      
    }) 
  })
  
  
  
  
  ############################# Differential expresssion #################################  
  observe({
    
    ResultsReactive()
    
  })
  
  
  ResultsReactive <- eventReactive(input$initResults,{
    withProgress(message = "Processing , please wait",{
      df2<-df_products()
      #remove the NA values from the counts dataframe  
      df2<-na.omit(df2) 
      
      #convert the counts dataframe to a data matrix
      countDataMatrix <- as.matrix(df2[ , -1])
      rownames(countDataMatrix) <- df2[ , 1]
      
      samples <- myValues$DF
      
      
      
      dds1 <- DESeqDataSetFromMatrix(countData=countDataMatrix,
                                     colData=samples,
                                     
                                     design=as.formula(input$designFormula))
      
      
      dds <- DESeq(dds1)
      myValues$dds<-dds
      
      contrast1<-input$condition1
      contrast2<-input$condition2
      
      myValues$countDataMatrix <- countDataMatrix
      myValues$sample_names<- colnames(countDataMatrix)
      
      #contrast column
      col<-input$factorNameInput
      #contrast column and levels
      contr <-  c(col, contrast1,contrast2)
      
      
      res1<- results(dds, contrast=contr, alpha= 0.05)
      myValues$res1<-res1
      
      
      #padj input
      padj.cutoff <- input$padj
      #lfc input
      lfc.cutoff <- input$log2F
      
      #convert the dds results to a dataframe and add the header name "gene_name"
      #to the gene ids column
      results <- res1 %>%
        data.frame() %>%
        rownames_to_column(var="gene_name") %>% 
        as_tibble()
      
      myValues$names<-rownames(res1)
      
      
      differ_expr <- results %>%
        filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)
      
      
      vsd <- varianceStabilizingTransformation(dds)
      myValues$vstMat <- assay(vsd)
      myValues$vsd <- vsd
      myValues$vsdColNames <- colnames(vsd)
      
      
      return(list('differ_expr'=differ_expr))
    })
  })
  
  
  
  
  #Read the DEGs file 
  df_products_upload3 <- reactive({
    inFile3 <- input$target_upload3
    
    if (is.null(inFile3))
      return(NULL)
    dataframe <- read.csv(inFile3$datapath, header = TRUE,sep = ',')
    #remove LOC pattern if you locate it
    if(sum(str_detect(dataframe$gene_name, 'LOC')) > 0){
      
      dataframe$gene_name<-sub("LOC","",as.character(dataframe$gene_name))}
    
    
    #Extra pattern removal(OPTIONAL)
    if(input$removePattern)
    {
      string<-input$pattern
      
      
      if(sum(str_detect(dataframe$gene_name, string)) > 0){
        
        dataframe$gene_name<-sub(string,"",as.character(dataframe$gene_name))}}
    
    
    
    return(dataframe)
    
    
  })
  
  
  output$sample_table<- DT::renderDataTable({
    df <- df_products_upload3()
    DT::datatable(df)
  })
  
  
  
  output$text <- renderPrint({
    
    
    if (!is.null(AnnotationReactive())) {
      
      return(myValues$annotation) 
      
    }
    
  })
  
  
  
  ######################## Annotation ############################# 
  observe({
    
    
    AnnotationReactive()
    
    
  })
  
  #search for the organism input in AnnotationHub data base
  AnnotationReactive <- eventReactive(input$search,{
    withProgress(message = "Searching , please wait",{
      organism <- input$organism
      
      
      hub <- AnnotationHub()
      annotation<-query(hub,c("OrgDb", organism))
      myValues$annotation<-annotation
      annotation
      
      
    })
    
  })
  
  #download the file that the user chose 
  observe({
    keytypesReactive()
  })
  keytypesReactive <- eventReactive(input$annotation,{
    withProgress(message = "Downloading , please wait",{
      # download the chosen file
      hub <- AnnotationHub()
      annot <- hub[[input$annot]]
      #update organism names for the transformation function 
      names<- genekitr::ensOrg_name
      latin_short_name<-names$latin_short_name
      updateSelectizeInput(session, "org", choices = latin_short_name, selected = latin_short_name[1],server=TRUE)
      
      myValues$annot<-annot
      annot
      
    })
    
    showNotification("Downloaded!", duration = 10, 
                     type = ( "message"))
  }) 
  
  
  
  
  
  
  
  
  
  ################################ Gene Set enrichment ############################
  myValues = reactiveValues()
  
  observe({
    
    
    gseGoReactive()
    
    
  })
  
  
  gseGoReactive <- eventReactive(input$initGo,{
    
    withProgress(message = "Processing , please wait",{
      
      isolate({
        
        
        
        validate(need(tryCatch({
          
          annotation <- myValues$annot
          df1<- df_products_upload3()
          
          #transform gene_names to ENTREZIDs
          f<-transId(df1$gene_name,org=input$org, transTo = "entrez")
          gene_name<-f$entrezid
          
          original_gene_list <- df1$log2FoldChange
          
          # name the vector
          names(original_gene_list) <- gene_name
          
          # omit any NA values 
          gene_list<-na.omit(original_gene_list)
          
          # sort the list in decreasing order (required for clusterProfiler)
          gene_list = sort(gene_list, decreasing = TRUE)
          
          myValues$gene_list<-gene_list
          
          
          setProgress(value = 0.3, detail = "Performing GSE analysis, please wait ...")
          
          go_gse <- gseGO(geneList=gene_list, 
                          ont = input$ontology, 
                          keyType = "ENTREZID",
                          #nPerm = input$nPerm, 
                          minGSSize = input$minGSSize, 
                          maxGSSize = input$maxGSSize, 
                          pvalueCutoff = input$pvalCuttoff, 
                          verbose = T, 
                          OrgDb = annotation, 
                          pAdjustMethod = input$pAdjustMethod)
          
          if(nrow(go_gse) < 1)
          {
            showNotification(id="warnNotify", "No gene can be mapped ...", type = "warning", duration = NULL)
            showNotification(id="warnNotify2", "Tune the parameters and try again.", type = "warning", duration = NULL)
            return(NULL)
          }
          
          updateNumericInput(session, "showCategory_dot", max = nrow(go_gse@result) , min = 0, value = ifelse(nrow(go_gse@result) > 0, 5,0))
          updateNumericInput(session, "showCategory_ridgeplot", max = nrow(go_gse@result) , min = 0, value = ifelse(nrow(go_gse@result) > 0, 5,0))
          updateNumericInput(session, "showCategory_cnet", max = nrow(go_gse@result) , min = 0, value = ifelse(nrow(go_gse@result) > 0, 5,0))
          updateNumericInput(session, "showCategory_gsea", max = nrow(go_gse@result) , min = 0, value = ifelse(nrow(go_gse@result) > 0, 5,0))
          updateNumericInput(session, "showCategory_heatplot", max = nrow(go_gse@result) , min = 0, value = ifelse(nrow(go_gse@result) > 0, 5,0))
          
        }, error = function(e) {
          myValues$status = paste("Error: ",e$message)
          
          showNotification(id="errorNotify", myValues$status, type = "error", duration = NULL)
          showNotification(id="errorNotify1", "Make sure that the gene_name column has the right format!", type = "error", duration = NULL)
          showNotification(id="errorNotify2", "Make sure that the needed  columns have the appropriate names! ", type = "error", duration = NULL)
          return(NULL)
        }
        
        ), 
        "Error merging files"))
        
        
      })
      
      
      return(list('go_gse'=go_gse))
      
      
    })
  }) 
  
  
  #Render the counts from the quantification process, as a data table 
  output$Table <- renderDataTable({
    featureCounts <- QuantifyReactive()
    
    if(!is.null(featureCounts)){
      df<-featureCounts$counts
      for ( col in 1:ncol(df)){
        colnames(df)[col] <-  sub("_trimmed", "", colnames(df)[col])
      }
      
      for ( col in 1:ncol(df)){
        colnames(df)[col] <-  sub("*.fastq.subread.SAM$", "", colnames(df)[col])
      }
      
      for ( col in 1:ncol(df)){
        colnames(df)[col] <-  sub("*.fastq.subread.BAM$", "", colnames(df)[col])
      }
      
      for ( col in 1:ncol(df)){
        colnames(df)[col] <-  sub("-", "_", colnames(df)[col])
      }
      
      
      for ( col in 1:ncol(df)){
        colnames(df)[col] <-  sub("R1", "", colnames(df)[col])
      }
      
      
      resultDF = df
    }
    
    
  })
  
  
  
  #Download the counts data table in .csv format
  output$downloadcountCSV <- downloadHandler(
    filename = function()  {paste0("counts",".csv")},
    content = function(file) {
      #Remove _trimmed string if you find it
      df=QuantifyReactive()$counts
      for ( col in 1:ncol(df)){
        colnames(df)[col] <-  sub("_trimmed", "", colnames(df)[col])
      }
      #Remove .fastq.subread.SAM if you find it
      for ( col in 1:ncol(df)){
        colnames(df)[col] <-  sub("*.fastq.subread.SAM$", "", colnames(df)[col])
      }
      #Remove .fastq.subread.BAM  string if you find it
      for ( col in 1:ncol(df)){
        colnames(df)[col] <-  sub("*.fastq.subread.BAM$", "", colnames(df)[col])
      }
      #Change - string to _ 
      for ( col in 1:ncol(df)){
        colnames(df)[col] <-  sub("-", "_", colnames(df)[col])
      }
      #Remove R1 string if you find it
      for ( col in 1:ncol(df)){
        colnames(df)[col] <-  sub("R1", "", colnames(df)[col])
      }
      
      write.csv(df, file, row.names=TRUE)}
  )
  
  output$Table2 <- renderDataTable({
    differ_expr<- ResultsReactive()$differ_expr
    
    if(!is.null(differ_expr)){
      
      resultDF =differ_expr
      
    }
    
  })
 
  
  #Download the differentially expressed genes table in .csv format 
  output$downloadDDSCSV2 <- downloadHandler(
    filename = function()  {paste0("dds",".csv")},
    content = function(file) {
      differ_expr<- ResultsReactive()$differ_expr
      
      
      write.csv( differ_expr, file, row.names=TRUE)}
    
  )
  
  #Render the enriched genes as a data table 
  output$gseGoTable <- renderDataTable({
    gseGo <- gseGoReactive()
    
    if(!is.null(gseGo)){
      resultDF = gseGo$go_gse@result
      
      DT::datatable(resultDF, options = list(scrollX = TRUE, columnDefs = list(list(visible=input$showAllColumns, targets= 10:12 )) ))
    }
    
  })
  
  #Download the enriched genes in .csv format 
  output$downloadgseGoCSV <- downloadHandler(
    filename = function()  {paste0("gsego",".csv")},
    content = function(file) {
      
      write.csv(gseGoReactive()$go_gse@result, file, row.names=TRUE)}
  )
  
  
  #Render Dotplot
  output$dotPlot = renderPlot({
    go_gse = gseGoReactive()$go_gse
    plot<-dotplot(go_gse, 
                  #number of enriched terms to display
                  showCategory = input$showCategory_dot, 
                  font.size = 8,
                  #separate result by 'category' variable
                  split=".sign") + facet_grid(.~.sign)
    
    
    myValues$plot<-plot
    print(plot)
    
  })
  
  #Download Dotplot
  output$downloadPlot <- downloadHandler(
    filename = function(){paste("dotplot", '.pdf', sep = '')},
    
    content = function(file){
      pdf(file, width = 8, height = 8)
      print(myValues$plot)
      dev.off()
    })
  
  
  
  #Render ridgePlot
  output$ridgePlot = renderPlot({
    go_gse = gseGoReactive()$go_gse
    ridge<- ridgeplot(go_gse, 
                      #number of enriched terms to display
                      showCategory = input$showCategory_ridge
    )
    
    myValues$ridge<-ridge
    print(ridge)
    
  })
  
  #Download ridgePlot
  output$downloadPlot2 <- downloadHandler(
    filename = function(){paste("ridge", '.pdf', sep = '')},
    
    content = function(file){
      pdf(file, width = 8, height = 8)
      print(myValues$ridge)
      dev.off()
    })
  
  
  
  
  #Render gseaplot
  output$gseaplot = renderPlot({
    go_gse = gseGoReactive()$go_gse
    gsea<- gseaplot(go_gse,
                    #"runningScore" and "position"
                    by = "all", 
                    title = go_gse$Description[input$geneSetId_gsea],
                    geneSetID = input$geneSetId_gsea)  
    
    
    myValues$gsea<-gsea
    print(gsea)
    
    
  })
  
  #Download gseaplot
  output$downloadPlot3 <- downloadHandler(
    filename = function(){paste("gsea", '.pdf', sep = '')},
    
    content = function(file){
      pdf(file, width = 8, height = 8)
      print(myValues$gsea)
      dev.off()
    })
  
  #Render cnetplot
  output$cnetplot = renderPlot({
    go_gse = gseGoReactive()$go_gse
    cnet<-cnetplot(go_gse,
                   categorySize="pvalue",
                   foldChange= myValues$gene_list,
                   #number of enriched terms to display
                   showCategory = input$showCategory_cne)
    
    myValues$cnet<-cnet
    print(cnet)
    
    
  })
  
  #Download cnetplot
  output$downloadPlot4 <- downloadHandler(
    filename = function(){paste("cnet", '.pdf', sep = '')},
    
    content = function(file){
      pdf(file, width = 8, height = 8)
      print(myValues$cnet)
      dev.off()
    })
  
  #Render heatplot
  output$heatplot = renderPlot({
    go_gse = gseGoReactive()$go_gse
    heat<-heatplot(go_gse, foldChange= myValues$gene_list)
    
    myValues$heat<-heat
    print(heat)
    
  })
  
  #Dowload heatplot 
  output$downloadPlot5 <- downloadHandler(
    filename = function(){paste("heatplot", '.pdf', sep = '')},
    
    content = function(file){
      pdf(file, width = 8, height = 8)
      print(myValues$heat)
      dev.off()
    })
  
  
  
  #Render maplot
  output$maplot2 = renderPlot({
    if(!is.null(ResultsReactive()))
    {
      de<- myValues$res1
      
      ma<-ggmaplot(de, 
                   #Accepted false discovery rate for considering genes as differentially expressed.
                   fdr = 0.05, 
                   #the fold change threshold.
                   #Only genes with a fold change >= fc and padj <= fdr are considered as significantly differentially expressed.
                   fc = 1, 
                   #points size
                   size = 0.6,
                   palette =c(input$col1,input$col2,input$col3),
                   genenames = as.vector(de$name),
                   legend = "top", top = 0,
                   #select.top.method = c( "fc"),
                   font.label = c("bold", 11),
                   font.legend = "bold",
                   font.main = "bold",
                   #ggplot2 theme name
                   ggtheme = ggplot2::theme_minimal())
      
      myValues$ma<-ma
      print(ma)
    }
  })
  
  #Download maplot
  output$downloadPlotma <- downloadHandler(
    filename = function(){paste("ma_plot", '.pdf', sep = '')},
    
    content = function(file){
      pdf(file, width = 8, height = 8)
      print(myValues$ma)
      dev.off()
    })
  
  
  
  #Render volcanoplot
  output$volcanoplot2 = renderPlot({
    if(!is.null(ResultsReactive()))
    {
      
      res1<- myValues$res1
      res2<- as.data.frame(res1)
      
      # add a column of NAs
      res2$diffexpressed <- "NO"
      # if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
      res2$diffexpressed[res2$log2FoldChange > 1 & res2$pvalue < 0.05] <- "UP"
      # if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
      res2$diffexpressed[res2$log2FoldChange < -1& res2$pvalue < 0.05] <- "DOWN"
      
      vol2<- ggplot(data=res2, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed)) +
        geom_point() + 
        theme_minimal() +
        #color choice and add lines
        scale_color_manual(values=c(input$col1,input$col3,input$col2)) +
        geom_vline(xintercept=c(-1, 1), col="black") +
        geom_hline(yintercept=-log10(0.05), col="black")
      
      
      myValues$vol2<-vol2
      print(vol2)
    }
  })
  
  
  #Download volcanoplot
  output$downloadPlotvol <- downloadHandler(
    filename = function(){paste("volcano_plot", '.pdf', sep = '')},
    
    content = function(file){
      pdf(file, width = 8, height = 8)
      print(myValues$vol2)
      dev.off()
    })
  
  
  
  #Update condition choices for pcaplot 
  observe({
    
    
    updateSelectInput(session,'group',
                      choices= colnames(myValues$DF),
                      selected = colnames(myValues$DF)[1])
    
  })
  
  #Render pcaplot 
  output$pcaplot = renderPlot({
    if(!is.null(ResultsReactive()))
    { 
      intgroups=input$group
      #pca function from DESeq2
      pca<- plotPCA(myValues$vsd, intgroup=intgroups)
      
      
      myValues$pca<-pca
      print(pca)
    }
  })
  
  #Download pcaplot
  output$downloadPlot7 <- downloadHandler(
    filename = function(){paste("pca", '.pdf', sep = '')},
    
    content = function(file){
      pdf(file, width = 8, height = 8)
      print(myValues$pca)
      dev.off()
    })
  
  
  
  #Update gene choices for boxplot
  observe({
    
    updateSelectizeInput(session,'sel_gene',
                         choices= myValues$names,
                         server=TRUE)
    
  })
  
  #Update condition choices for boxplot  
  observe({
    
    updateSelectInput(session,'boxplotX',
                      choices= colnames(myValues$DF),
                      selected = colnames(myValues$DF)[1])
    
    updateSelectInput(session,'boxplotFill',
                      choices= colnames(myValues$DF),
                      selected = colnames(myValues$DF)[1])
    
  })
  
  observe({
    geneExrReactive()
  })
  
  #
  geneExrReactive <- reactive({
    
    
    validate(need(length(input$sel_gene)>0,"Please select a gene."))
    #filter genes by log2
    filtered <- t(log2((counts(myValues$dds[input$sel_gene, ], normalized=TRUE, replaced=FALSE)+.5))) %>%
      merge(colData(myValues$dds), ., by="row.names") %>%
      gather(gene, expression, (ncol(.)-length(input$sel_gene)+1):ncol(.))
    
    
  })
  #Render boxPlot for gene expression 
  output$boxPlot <- renderPlot({
    if(!is.null(geneExrReactive()))
    {
      filtered = geneExrReactive()
      
      validate(need(length(input$boxplotX)>0,"Please select a group."))
      validate(need(length(input$boxplotFill)>0,"Please select a fill by group."))
      
      
      
      p <- ggplot(filtered, aes_string(input$boxplotX, "expression", fill=input$boxplotFill)) + 
        geom_boxplot() + scale_fill_manual(values=c(input$col1,input$col2))+facet_wrap(~gene, scales="free_y")
      
      
      
      myValues$p<-p
      print(p)
      
    }
    
    
  })
  
  #Download boxPlot for gene expression
  output$downloadPlot8 <- downloadHandler(
    filename = function(){paste("boxplot", '.pdf', sep = '')},
    
    content = function(file){
      pdf(file, width = 8, height = 8)
      print(myValues$p)
      dev.off()
    })
  
  
  
  
  #Update condition choices for heatmap
  observe({
    
    updateSelectInput(session,'group1',
                      choices= colnames(myValues$DF),
                      selected = colnames(myValues$DF)[1])
    
    updateSelectInput(session,'group2',
                      choices= colnames(myValues$DF),
                      selected = colnames(myValues$DF)[1])
    
  })
  
  
  observe({
    heatmapReactive()
  })
  
  heatmapReactive <- reactive({
    if(!is.null(ResultsReactive()))
      
      #for 1 condition
      if (input$m == 1) 
      {
        isolate({
          
          df <-as.data.frame(colData(myValues$dds)[,(input$group1)])
          
          row.names(df) <- colnames(myValues$vsd)
          colnames(df)[1] <- input$group1
          #select the top 30 differentially expressed genes 
          topVarGenes <- head(order(-rowVars(assay(myValues$vsd))),30)
          myValues$df<-df
          mat <- assay(myValues$vsd)[ topVarGenes, ]
          mat <- mat - rowMeans(mat)
          
          return(mat)
        })
      }
    #for 2 conditions 
    if (input$m == 2) { 
      isolate({
        
        df <-as.data.frame(colData(myValues$dds)[,c(input$group1,input$group2)])
        row.names(df) <- colnames(myValues$vsd)
        #select the top 30 differentially expressed genes 
        topVarGenes <- head(order(-rowVars(assay(myValues$vsd))),30)
        myValues$df<-df
        
        mat <- assay(myValues$vsd)[ topVarGenes, ]
        mat <- mat - rowMeans(mat)
        
        return(mat)
      })
    }
    
    
    
  })  
  
  
  #Render heatmap
  output$heatmap <- renderPlot({
    if(!is.null(heatmapReactive()))
    {
      mat=heatmapReactive()
      
      
      heatmap<- pheatmap(mat,
                         annotation_col=myValues$df,
                         scale = "row", 
                         cluster_rows = T,
                         angle_col = 90,
                         
                         
      )
      
      myValues$heatmap<-heatmap
      print(heatmap)
      
    }
  })
  
  
  #Download heatmap
  output$downloadPlot6 <- downloadHandler(
    filename = function(){paste("heatmap", '.pdf', sep = '')},
    
    content = function(file){
      pdf(file, width = 8, height = 8)
      print(myValues$heatmap)
      dev.off()
    })
  
  
  #Render boxplot for counts
  output$boxplot2 <- renderPlot({
    if(!is.null(myValues$dds))
    {
      
      min_nonzero=1
      #filtered counts by log2
      b<-boxplot(log2(myValues$countDataMatrix+min_nonzero),
                 col= brewer.pal(n = ncol(myValues$countDataMatrix),
                                 name = input$Palette),
                 names=myValues$sample_names, 
                 las=2, 
                 ylab="log2(counts)",
                 main="Distribution of counts for all samples")
      myValues$b<-b
      print(b)
      
    }
  })
  
  
  
  
  
  
}
shinyApp(ui=ui, server=server)



