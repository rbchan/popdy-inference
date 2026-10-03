library(shiny)
library(rsconnect)

## To deploy the app, follow instructions here:
## https://docs.rstudio.com/shinyapps.io/getting-started.html#deploying-applications


if(1==2) {
    ## If getwd() is same as script dir
    runApp()
    deployApp(appName="scr-cap-prob")

    ## If getwd() is up one level
    runApp("shiny")
    deployApp("shiny", appName="scr-cap-prob")



## Move up one dir level
rsconnect::migrateToConnectCloud(
  appPath = "shiny",
  contentId = "01a10397-4d53-ad4c-d21e-f48c0e57e7df"
)


}
