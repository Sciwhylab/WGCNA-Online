# Shiny app docker file
# Based on example at
# https://blog.sellorm.com/2021/04/25/shiny-app-in-docker/

# get shiny server and R from the rocker project
FROM rocker/shiny

# system libraries
# Try to only install system libraries you actually need
# Package Manager is a good resource to help discover system deps
RUN apt-get update && apt-get upgrade -y

# install R packages required 
# Change the packages list to suit your needs
RUN R -e 'install.packages(c(\
              "shiny", \
              "bslib", \
              "BiocManager", \
              "DT"\
            )\
          )'

# RUN R -e 'BiocManager::install(c(\
#               "impute", \
#               "preprocessCore", \
#               "GO.db", \
#               "AnnotationDbi"\
#               )\
#             )'

RUN R -e 'BiocManager::install(c(\
              "WGCNA", \
              "impute", \
              "flashClust" \
            )\
          )'

RUN R -e 'BiocManager::valid()'

# copy the app directory into the image
COPY ./* /srv/shiny-server/RefEx_Expression/

# run app
CMD ["/usr/bin/shiny-server"]
