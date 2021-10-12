pageWithSidebar(
    titlePanel('Distance sampling (half-normal detection function)'),
    sidebarPanel(
        selectInput('transect', 'Survey type', c("Line transect", "Point transect")),
        numericInput('sigma', 'Scale parameter', 20, min=0, max=1000, step=5),
        numericInput('distanceMax', 'Maximum distance', 100, min=0, max=500, step=10),
    ),
    mainPanel(
        fluidRow(
            column(6, plotOutput('gx')),
            column(6, plotOutput('px'))
        ),
        fluidRow(
            column(6, plotOutput('gxpx')),
            column(6,
                h4("Average detection prob\n"),
                tableOutput('pbar')
                )
        )
    )
)
