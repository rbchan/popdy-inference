# popdy-inference
A course on statistical inference for models of population dynamics. The course is aimed at wildlife ecology graduate students at UGA. 

## Course objectives
Students will learn how to fit statistical models to field data and draw inferences on the processes governing spatial and temporal variation in occupancy and abundance. Focus will be on occupancy data, count data, distance sampling data, mark-recapture data, and telemetry data. Hierarchical statistical models will be used for inference.

## Prerequisites
Students should have taken a good introductory course on population dynamics, and they should be proficient in [R](https://www.r-project.org/). 

## Structure
Class meetings will consist of a mixture of lectures and computer exercises.

## Assignments
Students will analyze their own data, write up their results, and present their findings at the end of the semester. Two rounds of peer review will be involved. Weekly computer problems will also be assigned. 

## Building the PDF slides
Lecture slides created using [LaTeX](https://www.latex-project.org/), [Beamer](https://en.wikipedia.org/wiki/Beamer_(LaTeX)), and [knitr](https://yihui.org/knitr/). To build the PDFs, open [R]((https://www.r-project.org/) and navigate to a lecture directory, such as `lectures/Nmix-binII` and issue the following commands:

```
library(knitr)
knit("lecture-Nmix-binomial-II.Rnw")
```
This will produce a standard Beamer `.tex` file that you can compile as you would any other Latex file. 


## History
This course is based loosely off of "Estimation of Fish and Wildlife Population Parameters (WILD 8390)", a course taught by Dr. Mike Conroy for approximately 30 years. I had the great fortune of co-instructing the course with Mike for two years. 
