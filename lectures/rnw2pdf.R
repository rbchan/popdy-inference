rnw2pdf <- function(filestub, env=globalenv(), tangle=FALSE,
                    open=TRUE, clean=TRUE, ...) {

    ## Load knitr or install if need be
    reqval <- require(knitr)
    if(!reqval) {
        install.packages("knitr")
        library(knitr)
    }

    ## Create the .tex file
    rnw.file <- paste(filestub, ".Rnw", sep="")
    tex.file <- knit(rnw.file, envir=env, ...)

    ## Create the .R file
    if(tangle)
        knit(rnw.file, tangle=TRUE)

    ## Create the PDF
    ## Make sure to install texinfo and texlive on linux
    tools::texi2dvi(tex.file, pdf=TRUE, clean=clean)

    ## Open the PDF
    if(open) {
        os <- .Platform$OS.type
        if(os=="unix") {
##            open.pdf <- paste0("okular ", filestub, ".pdf") ## Must install okular
            open.pdf <- paste0("gopen ", filestub, ".pdf") ## Must install gopen
##            open.pdf <- paste0("xdg-open ", filestub, ".pdf")
        } else 
            open.pdf <- paste0("open ", filestub, ".pdf") ## Requires Rtools
        sysval <- system(open.pdf, wait=FALSE, ignore.stdout=TRUE, ignore.stderr=TRUE)
        if(sysval != 0) {
            warning("\n\nThe PDF should be in your working directory")
        }
    }

}
