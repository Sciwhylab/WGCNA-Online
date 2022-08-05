# WGCNA Online

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
