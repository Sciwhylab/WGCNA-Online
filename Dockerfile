# Use an official R runtime as a parent image
# https://rocker-project.org/images/versioned/shiny.html
FROM rocker/shiny:latest

# Install any needed packages
RUN R -e "install.packages(c(\"shiny\", \"shinyWidgets\", \"stats\", \"bslib\", \"BiocManager\", \"DT\", \"here\", \"readr\", \"readxl\"))"
RUN R -e "BiocManager::install(c(\"WGCNA\", \"impute\", \"flashClust\"))"

# Set the working directory in the container to /app
WORKDIR /app

# Copy the current directory contents into the container at /app
COPY . /app

# Run app.R when the container launches
CMD ["R", "-e", "shiny::runApp('/app', host = '0.0.0.0', port = 3838)"]
