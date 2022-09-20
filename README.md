# WGCNA Online

An online application that lets users perform weighted gene co-expression network analysis on their data –  
a no-code solution that helps you quickly get started and explore.

This project aims to create a GUI interface to the [R WGCNA package](https://www.rdocumentation.org/packages/WGCNA/versions/1.71).

## Live version

A [hosted version](https://kitswas.shinyapps.io/WGCNA-Online/) is available.  
Note that it is _resource-limited._

## Getting started

### Prerequisites

[R](https://cloud.r-project.org/) and [git](https://git-scm.com/downloads) must be installed on your system.

### Clone the repository

```shell
git clone https://github.com/Sciwhylab/WGCNA-Online.git
cd WGCNA-Online
```

### Install required packages

In R, execute the following commands:

```r
install.packages(c("shiny", "bslib", "BiocManager", "DT"))
```

```r
BiocManager::install(c("WGCNA", "impute", "flashClust"))
```

If you are having trouble with the above steps, refer to this [answer on StackOverflow](https://stackoverflow.com/a/50364335/8659747).

### Run the Shiny App

First clone the repository.  
Then, go to the project root and execute the following in R.

```r
shiny::runApp(, host = "0.0.0.0")
```

## Troubleshooting

If you are facing package installation problems while deploying the app to [shinyapps.io](https://www.shinyapps.io/), execute the following in R.

```r
options(repos = c(getOption("repos"),  BiocManager::repositories()))
```
