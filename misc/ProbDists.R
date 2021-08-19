#need to add log link functions

library(shiny)
library(extraDistr)
ui<-navbarPage("Stats are Fun!",
    tabPanel("Probability Distributions",           
    fluidPage(
  titlePanel("Probability Distribution", windowTitle = "Probability Dists"),
  fluidRow(
    column(2,
           selectInput("dist","Choose the probability distribution: ", choices =
                         list(Continuous = list("Normal", "Half-Normal", "Beta", "Gamma", "Uniform"),
                              Discrete = list("Binomial", "Poisson")))
    ),
    column(2,
           conditionalPanel(
             condition = "input.dist == 'Normal'",
             numericInput("mean","Mean: ", value = 0),
             numericInput("dev", "Standard deviation: ", value = 1,min=0)
           ),
           conditionalPanel(
             condition = "input.dist == 'Half-Normal'",
             numericInput("sigma", "Sigma: ", value = 15, min = 0)
           ),
           conditionalPanel(
             condition = "input.dist == 'Gamma'",
             numericInput("shape", "Shape (alpha): ", value = 1, min = 0),
             numericInput("rate", "Rate (beta): ", value = 1, min = 0)
             
           ),
           conditionalPanel(
             condition = "input.dist == 'Uniform'",
             numericInput("a", "a", value = 0),
             numericInput("b", "b", value = 1),
             
           ),
           conditionalPanel(
             condition = "input.dist == 'Beta'",
             numericInput("alpha","Alpha: ", value = 1,min=0),
             numericInput("beta", "Beta: ", value = 1,min=0)
           ),
           conditionalPanel(
             condition = "input.dist == 'Poisson'",
             numericInput("lambda","Lambda: ", value = 9,min = 0)
           ),
           conditionalPanel(
             condition = "input.dist == 'Binomial'",
             numericInput("n","n: ", value = 10,step=1),
             numericInput("p", "p: ", value = 0.5, min = 0, max = 1)
           ),
    ),
    
  ),
  fluidRow(
    column(6,
           plotOutput("plot")
    ),
    column(6,
           helpText("You can easily graph these distributions in R!"),
           helpText("First make a vector of points to plot (ex: -30,30), then use the built in R fuctions to calculate the density!"),
           conditionalPanel(
             condition = "input.dist == 'Normal'",
                helpText(code("points <- seq(from=-30,to=30,length.out=1000)")),
                helpText(code("plot(points, dnorm(points,mean,sd), type = 'l')"))),
           conditionalPanel(
             condition = "input.dist == 'Beta'",
             helpText(code("points <- seq(from=-30,to=30,length.out=1000)")),
                helpText(code("plot(points, dbeta(points,alpha,beta), type = 'l'"))),
           conditionalPanel(
             condition = "input.dist == 'Gamma'",
             helpText(code("points <- seq(from=-30,to=30,length.out=1000)")),
                helpText(code("plot(points, dgamma(points, shape, rate), type = 'l')"))),
           conditionalPanel(
             condition = "input.dist == 'Uniform'",
             helpText(code("points <- seq(from=-30,to=30,length.out=1000)")),
           helpText(code("plot(points, dunif(points, min, max), type = 'l')"))),
           conditionalPanel(
             condition = "input.dist == 'Half-Normal'",
             helpText("For the Half-Normal distribution you also need to load a package called extraDistr"),
             helpText(code("library(extraDistr)")),
             helpText(code("points <- seq(from=-30,to=30,length.out=1000)")),
             helpText(code("plot(points, dhnorm(points,sigma), type = 'l')"))),
           conditionalPanel(
             condition = "input.dist == 'Binomial'",
             helpText("For discrete distributions, you'll want to add", code("type = 'h'"), "to your plotting function"),
             helpText(code("points <- seq(from=-30,to=30,length.out=1000)")),
             helpText(code("plot(points, dbinom(points,n,p), type = 'h')"))),
           conditionalPanel(
             condition = "input.dist == 'Poisson'",
             helpText("For discrete distributions, you'll want to add", code("type = 'h'"), "to your plotting function"),
             helpText(code("points <- seq(from=-30,to=30,length.out=1000)")),
             helpText(code("plot(points, dpois(points,lambda), type = 'h')")))
           
          
           )
))),
tabPanel("Logistic Regression",
  fluidPage(
    titlePanel("Logistic Regression", windowTitle = "Logistic Regression"),
      fluidRow(
         column(2, 
            numericInput("int", "Intercept", value = 5, step = .1),
            numericInput("slope", "Slope", value = -1, step = .1)),
         column(3,
                plotOutput("plot_lin")
         ),
         column(3,
                plotOutput("plot_logist")
         ),
         column(3,
                plotOutput("plot_odds")
         )
         ),
    fluidRow(
      column(2, 
             numericInput("x", "Predict value when X = ", value = 5, step = .1)),
      column(3, 
             textOutput("pred_lin") ),
      column(3,
             textOutput("pred_logistic")),
      column(3,
             textOutput("pred_odds"))),
    fluidRow(
      column(3),
      column(4,
             helpText("Predicted Value =", code("int+X*slope"), "on the linear scale"),
             helpText("This is also ln(p/1-p), though of course ")),
      column(4, 
             textOutput("code_logistic"),
             helpText("Predicted Value =", code("plogis(int+X*slope)"), "on the probability scale"),
             helpText("Remember that plogis is just the inverse logit function, so you can alternatively write", code("1/(1+exp(-(int+X*slope))")))
    )# fluid row end
    ) #fluid page end
  )# tab panel end
) #end ui

server<-function(input, output, session){
#Panel 1
  dist <-reactive(input$dist)
  mean <-reactive(input$mean)
  sd <-reactive(input$dev)
  lambda <-reactive(input$lambda)
  n <-reactive(input$n)
  p <-reactive(input$p)
  alpha <-reactive(input$alpha)
  beta <-reactive(input$beta)
  sigma <-reactive(input$sigma)
  a <-reactive(input$a)
  b <- reactive(input$b)
  shape <-reactive(input$shape)
  rate <-reactive(input$rate)
  x <- reactive(input$x)
  
  output$plot<-renderPlot({
    
    lim1<-switch(dist(),"Normal"=-30,"Beta"=0,"Poisson"=0,"Binomial"=0,"Half-Normal"=0, "Uniform" = a(), "Gamma" = 0)
    lim2<-switch(dist(),"Normal"=30,"Beta"=1,"Poisson"=2*lambda(),"Binomial"=n(),"Half-Normal"=50, "Uniform" = b(), "Gamma" = 20)
    points<-switch(dist(),"Normal"=seq(from=lim1,to=lim2,length.out=1000),"Beta"=seq(from=lim1,to=lim2,length.out=1000), "Gamma" = seq(from=lim1,to=lim2,length.out=1000), "Uniform" = seq(from=lim1,to=lim2,length.out=1000),
                   "Poisson"=lim1:lim2,"Binomial"=lim1:lim2, "Half-Normal"=seq(from=lim1,to=lim2,length.out = 1000))
    Density<-switch(dist(),"Normal"=dnorm(points,mean(),sd()),
                    "Beta"=dbeta(points,alpha(),beta()),
                    "Poisson"=dpois(points,lambda()),
                    "Binomial"=dbinom(points,n(),p()),
                    "Half-Normal"=dhnorm(points,sigma()),
                    "Gamma" = dgamma(points, shape(), rate()),
                     "Uniform" = dunif(points, a(), b()))
    cx <- switch(dist(), "Normal" = .1, "Beta" = .1, "Poisson" = 1, "Binomial" = 1, "Half-Normal" = .1,"Gamma" = .1,"Uniform" = .1)
    tt <- switch(dist(), "Normal" = 'l', "Beta" = 'l', "Poisson" = 'h', "Binomial" = 'h', "Half-Normal" = 'l',"Uniform" = 'l',"Gamma" = 'l')
    mn<-switch(dist(),"Normal"=mean(),"Beta"=alpha()/(alpha()+beta()),"Poisson"=lambda(),"Binomial" = n()*p(),"Half-Normal" = sigma()*sqrt(2)/sqrt(3.14159), "Uniform" = (a()+b())/2, "Gamma" = shape()/rate())
    yl<<-c(0,4*max(Density)/3)
    plot(points,Density,pch = 19, col = "cyan3", type = tt, cex = cx, lwd=2.5, xlab="Values",main=dist(),ylim=yl)
    abline(v = mn, col="grey",lwd = 2, lty = 2)})#plot panel 1

#Panel 2
int  <- reactive(input$int)
slope <- reactive(input$slope)

output$plot_lin<-renderPlot({
  plot(seq(0,30, by = .1), int()+seq(0,30,by=.1)*slope(), type = "l", xlim = c(0,30), ylim = c(-100,100), ylab = "Predicted Value", xlab = "Covariate Value")
  points(x(), int()+x()*slope(), cex = 1.2, pch = 19, col = "red")
  
}) # linear regression 

output$plot_logist<-renderPlot({
  plot(seq(0,30, by = .1), plogis(int()+seq(0,30,by=.1)*slope()), type = "l", xlim = c(0,30), ylim = c(0,1), ylab = "Predicted Probability", xlab = "Covariate Value")
  points(x(), plogis(int()+x()*slope()), cex = 1.2, pch = 19, col = "red")
  
})# logistic regression 
output$plot_odds<-renderPlot({
  plot(seq(0,30, by = .1), plogis(int()+seq(0,30,by=.1)*slope())/(1-plogis(int()+seq(0,30,by=.1)*slope())), type = "l", xlim = c(0,30), ylim = c(0,1000), ylab = "Odds Ratio p/(1-p)", xlab = "Covariate Value")
  points(x(), plogis(int()+x()*slope())/(1-plogis(int()+x()*slope())), cex = 1.2, pch = 19, col = "red")
  
}) #plot of the odds; kind of useless haha





output$pred_lin <- renderText(paste("Predicted Value =", round(int()+x()*slope(), digits = 3) ))
output$pred_logistic <- renderText(paste("Predicted Probability =", round(plogis(int()+x()*slope()), digits = 3)))
output$pred_odds <- renderText(paste("Predicted Odds Ratio =", round(plogis(int()+x()*slope())/(1-plogis(int()+x()*slope())), digits = 3), "Note: Not very useful"))


}

shinyApp(ui, server)
