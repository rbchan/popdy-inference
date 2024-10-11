function(input, output, session) {

    output$hpp <- renderPlot({
        ED <- input$ED
        input$simButton
        N <- rpois(n=1, lambda=ED)
        s <- cbind(runif(N), runif(N))
        plot(s, pch=16, col="blue", asp=1, xlim=c(0,1), ylim=c(0,1),
##             main="Homogeneous Poisson point process",
             main=paste("N =", N),
             xlab="x coordinate", ylab="y coordinate")
    })
    
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

    output$spider <- renderPlot({
      ED <- input$EDFull
      input$simButtonFull
      N <- rpois(n=1, lambda=ED)
      s <- cbind(runif(N, 0, 1),
                 runif(N, 0, 1))
      g0 <- input$g0Full
      sig <- input$sigmaFull
      x <- cbind(rep(seq(0.15, 0.85, by=0.1), each=8),
                 rep(seq(0.15, 0.85, by=0.1), times=8))  ## Trap locations
      J <- nrow(x)                 ## nTraps
      dist.sx <- matrix(NA, N, J)  
      for(i in 1:N) {
        dist.sx[i,] <- sqrt((s[i,1]-x[,1])^2 + (s[i,2]-x[,2])^2)
      }
      p <- g0*exp(-dist.sx/(2*sig*sig))
      K <- 5                          # nOccasions
      y.all <- array(NA, c(N, J, K))
      for(i in 1:N) {
        for(j in 1:J) {
          y.all[i,j,] <- rbinom(K, 1, prob=p[i,j])
        }
      }
      captured <- rowSums(y.all)>0
      y <- y.all[captured,,]
      
      plot(0, xlim=c(0,1), ylim=c(0,1), asp=1, xlab="", ylab="",
           main="Activity centers, traps, and capture locs")
      s.cap <- s[captured,]
      for(i in 1:nrow(y)) {
        traps.i <- which(rowSums(y[i,,])>0)
        for(j in 1:length(traps.i)) {
          segments(s.cap[i,1], s.cap[i,2],
                   x[traps.i[j],1], x[traps.i[j],2], col=gray(0.3))
        }
      }
      points(s[captured,], pch=16, col="blue") ## Activity center locations
      points(s[!captured,], pch=1, col="blue") ## Activity center locations
      points(x, pch=3)                         ## Trap locations
      
    })
}
