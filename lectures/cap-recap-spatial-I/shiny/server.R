function(input, output, session) {

    output$gdist <- renderPlot({
        g0 <- input$g0
        sig <- input$sigma
        plot(function(x) g0*exp(-x^2 / (2*sig^2)), from=0, to=input$buffer,
             ylim=c(0, 1), col="blue", lwd=2,
             xlab="Distance (x)", ylab="g(x)", main="Capture probability")
    })

    output$gsp <- renderPlot({
        g0 <- input$g0
        sig <- input$sigma
        buffer <- input$buffer
        grid0 <- seq(from=0, to=buffer, by=input$resolution)
        grid <- cbind(rep(grid0, each=length(grid0)),
                      rep(grid0, times=length(grid0)))
        x <- matrix(c(buffer/2, buffer/2), nrow=1)
        distSq <- (grid[,1]-x[1])^2 + (grid[,2]-x[2])^2
        p <- g0*exp(-distSq/(2*sig*sig))
        library(lattice)
        levelplot(p ~ grid[,1]+grid[,2], aspect="iso",
                  panel=function(...) {
                      panel.levelplot(...)
                      panel.xyplot(x[,1], x[,2], pch=3, col="black")
                  }, 
                  xlab="x coordinate", ylab="y coordinate", main="Capture probability")
    })

    
}
