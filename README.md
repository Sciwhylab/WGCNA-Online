# WGCNA Online

An online application that lets users perform weighted gene co-expression network analysis on their data –  
a no-code solution that helps you quickly get started and explore.

This project aims to create a GUI interface to the [R WGCNA package](https://www.rdocumentation.org/packages/WGCNA/versions/1.71).

## Live version

A [hosted version](https://kitswas.shinyapps.io/WGCNA-Online/) is available on shinyapps.io. Note that it is _resource-limited._  
The app sleeps to save resources when idle, leading to longer startup times. Please be patient.

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
install.packages(c("shiny", "shinyWidgets", "stats", "bslib", "BiocManager", "DT", "here", "readr", "readxl"))
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

## Using the Dockerfile

If you have [Docker](https://www.docker.com/products/docker-desktop) installed, you can build and run the app using the provided Dockerfile.

You can build the image using the following command:

```bash
docker build -t wgcna-online .
```

And run the container using: (replace `8080` with the port you want to use)

```bash
docker run --name "wgcna" --rm -p 8080:3838 wgcna-online
```

Visit [http://localhost:8080](http://localhost:8080) to access the app.  
If you are using a server, replace `localhost` with the server's IP address.

## Troubleshooting

If you are facing package installation problems while deploying the app to [shinyapps.io](https://www.shinyapps.io/), execute the following in R.

```r
options(repos = c(getOption("repos"),  BiocManager::repositories()))
```
