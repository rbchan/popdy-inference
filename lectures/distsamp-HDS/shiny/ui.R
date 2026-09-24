## pageWithSidebar(
fluidPage(
    titlePanel('Distance sampling (half-normal detection function)'),
    sidebarPanel(
        selectInput('transect', 'Survey type', c("Line transect", "Point transect")),
        numericInput('sigma', 'Scale parameter', 20, min=0, max=1000, step=5),
        numericInput('distanceMax', 'Maximum distance', 100, min=0, max=500, step=10),
    width=3),
    mainPanel(
        fluidRow(
            column(4, plotOutput('gx', inline=TRUE)),
            column(4, plotOutput('px', inline=TRUE)),
            column(4, plotOutput('gxpx', inline=TRUE))
        ),
        fluidRow(
            column(4),
            column(4, 
                   h4("Average detection prob\n"),
                   tableOutput('pbar'),
            column(4),
            )
        )
    )
)
