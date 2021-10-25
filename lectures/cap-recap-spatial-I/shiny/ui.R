## pageWithSidebar(
fluidPage(
    titlePanel('SCR encounter function (Gaussian)'),
    tabsetPanel(
        ##    tabPanel("Plot", plotOutput("plot")),
        tabPanel("Summary", verbatimTextOutput("summary")),
        ##    tabPanel("Table", tableOutput("table")),
        tabPanel("Table", verbatimTextOutput("table"))
    ),
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
)
