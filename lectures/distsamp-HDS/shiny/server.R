function(input, output, session) {

    output$gx <- renderPlot({
        sig <- input$sigma
        plot(function(x) exp(-x^2 / (2*sig^2)), from=0, to=input$distanceMax,
             ylim=c(0, 1), xlab="Distance (x)", ylab="g(x)", main="Detection probability")
    })

    output$px <- renderPlot({
        plot(function(x) 1/rep(input$distanceMax, length(x)),
             from=0, to=input$distanceMax, ylim=c(0, 2/input$distanceMax),
             xlab="Distance (x)", ylab="p(x)", main="p(animal occurs at distance x)")
    })

    output$gxpx <- renderPlot({
        sig <- input$sigma
        plot(function(x) exp(-x^2 / (2*sig^2))/input$distanceMax,
             from=0, to=input$distanceMax, ylim=c(0, 1/input$distanceMax),
             xlab="Distance (x)", ylab="g(x)p(x)",
             main="p(animal occurs and is detected at distance x)")
    })

    output$pbar <- renderTable({
        sig <- input$sigma
        g <- function(x) exp(-x^2/(2*sig^2))    ## g(x)
        pdf <- function(x) 1/input$distanceMax  ## p(x)
        gp <- function(x) g(x)*pdf(x)
        pbar <- integrate(gp, lower=0, upper=input$distanceMax)$value
        pbarOut <- rbind("pbar" = pbar)
        return(pbarOut)
    }, rownames=TRUE, colnames=FALSE)
    
    
}
