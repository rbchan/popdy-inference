## pageWithSidebar(
navbarPage('SCR',
    ## titlePanel('SCR state and observation models'),
    ## tabsetPanel(
    ##     ##    tabPanel("Plot", plotOutput("plot")),
    ##     tabPanel("Summary", verbatimTextOutput("summary")),
    ##     ##    tabPanel("Table", tableOutput("table")),
    ##     tabPanel("Table", verbatimTextOutput("table"))
    ## ),
##    tabsetPanel(
  tabPanel("State model", ##verbatimTextOutput("State")),
           titlePanel("Homogeneous Poisson Point Process"),
           
           sidebarLayout(
             sidebarPanel(
               numericInput('ED', 'Expected value of density', 10, min=0, max=1000, step=1),
               actionButton('simButton', 'Simulate', class = "btn-success"),
               ),
             mainPanel(
               fluidRow(
                 column(7, plotOutput('hpp'))
               )
             )
           )),
  tabPanel("Observation model", 
           sidebarLayout(
             sidebarPanel(
               numericInput('g0', 'Baseline capture probability', 0.5, min=0, max=1, step=0.02),
               numericInput('sigma', 'Scale parameter', 20, min=0, max=1000, step=5),
               numericInput('buffer', 'Buffer (square)', 100, min=0, max=500, step=10),
               numericInput('resolution', 'Spatial resolution', 1, min=0, max=10, step=0.1),
               ),
             mainPanel(
               fluidRow(
                 column(5, plotOutput('gdist')),
                 column(7, plotOutput('gsp'))
               )
             )
           )
           ),
  ##    )
  tabPanel("Full model", 
           sidebarLayout(
             sidebarPanel(
               numericInput('EDFull', 'Expected value of density', 50, min=5, max=1000, step=1),
               numericInput('g0Full', 'Baseline capture probability', 0.1, min=0, max=1, step=0.01),
               numericInput('sigmaFull', 'Scale parameter', 0.15, min=0.01, max=0.5, step=0.01),
               actionButton('simButtonFull', 'Simulate', class = "btn-success"),
               ),
             mainPanel(
               fluidRow(
                 column(7, plotOutput('spider'))
               )
               ## plotOutput('spider')
             )
           )
           )
  )
