library(shiny)
library(rsconnect)

## To deploy the app, follow instructions here:
## https://docs.rstudio.com/shinyapps.io/getting-started.html#deploying-applications


if(1==2) {
    ## If getwd() is same as script dir
    runApp()
    deployApp(appName="distance-sampling")

    ## If getwd() is up one level
    runApp("shiny")
    deployApp(appDir="shiny", appName="distance-sampling")
}


if(1==3) {

## Move up one dir level
rsconnect::migrateToConnectCloud(
  appPath = "shiny",
  contentId = "01a0d342-62e7-2e67-b9bc-3d69064dd322"
)

## rsconnect::migrateToConnectCloud(
##   appPath = "shiny",
##   contentId = "abc123",
##   # Pass only the argument(s) needed to single out the record;
##   # often appName alone is enough.
##   appName = "my-dashboard",
##   account = "my-shinyapps-account",
##   server = "shinyapps.io"
## )    

}
