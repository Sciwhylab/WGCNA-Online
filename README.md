RefEx_expression
================

# Getting started

## Prerequisites

[R](https://cloud.r-project.org/) and [git](https://git-scm.com/downloads) must be installed on your system.

## Install required packages

```
install.packages(c("shiny", "bslib", "BiocManager", "DT"))
```

```
BiocManager::install(c("WGCNA", "impute", "flashClust"))
```

If you are having trouble with the above steps, refer to this [answer on StackOverflow](https://stackoverflow.com/a/50364335/8659747).

## Run the Shiny App

First clone the repository.  
Then, go to the project root and execute
```
shiny::runApp(, host = "0.0.0.0")
```
