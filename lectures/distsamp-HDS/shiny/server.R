function(input, output, session) {

    output$gx <- renderPlot({
        sig <- input$sigma
        plot(function(x) exp(-x^2 / (2*sig^2)), from=0, to=input$distanceMax,
             ylim=c(0, 1), col="blue", lwd=2,
             xlab="Distance (x)", ylab="g(x)",
             main="Detection probability", cex.main=1)
    }, width=252, height=252)

    output$px <- renderPlot({
        if(input$transect=="Line transect") {
            plot(function(x) 1/rep(input$distanceMax, length(x)),
                 from=0, to=input$distanceMax, ylim=c(0, 2/input$distanceMax),
                 col="orange", lwd=2,
                 xlab="Distance (x)", ylab="p(x)",
                 main="p(occurs at distance x)", cex.main=1)
        } else {
            plot(function(x) 2*x/input$distanceMax^2,
                 from=0, to=input$distanceMax, ylim=c(0, 2/input$distanceMax),
                 col="orange", lwd=2,
                 xlab="Distance (x)", ylab="p(x)",
                 main="p(occurs at distance x)", cex.main=1)
        }
    }, width=252, height=252)

    output$gxpx <- renderPlot({
        sig <- input$sigma
        if(input$transect=="Line transect") {
            plot(function(x) exp(-x^2 / (2*sig^2))/input$distanceMax,
                 from=0, to=input$distanceMax, ylim=c(0, 1/input$distanceMax),
                 col="seagreen", lwd=2,
                 xlab="Distance (x)", ylab="g(x)p(x)",
                 main="p(occurs and detected at distance x)", cex.main=1)
        } else {
            plot(function(x) exp(-x^2 / (2*sig^2))*(2*x/input$distanceMax^2),
                 from=0, to=input$distanceMax, ylim=c(0, 2/input$distanceMax),
                 col="seagreen", lwd=2,
                 xlab="Distance (x)", ylab="g(x)p(x)",
                 main="p(occurs and detected at dist x)", cex.main=1)
        }
    }, width=252, height=252)

    output$pbar <- renderTable({
        sig <- input$sigma
        g <- function(x) exp(-x^2/(2*sig^2))            ## g(x)
        if(input$transect=="Line transect") 
            pdf <- function(x) 1/input$distanceMax      ## p(x) 
        else 
            pdf <- function(x) 2*x/input$distanceMax^2  ## p(x)
        gp <- function(x) g(x)*pdf(x)
        pbar <- integrate(gp, lower=0, upper=input$distanceMax)$value
        pbarOut <- rbind("pbar" = pbar)
        return(pbarOut)
    }, rownames=TRUE, colnames=FALSE, bordered=TRUE)
    
    
}
