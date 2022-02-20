# DOCKERFILE for ABUNDANCER
# github.com/davidsbutcher/abundancer
# docker build --tag davidsbutcher/abundancer:test01 --no-cache .
# docker run -p 3838:3838 davidsbutcher/abundancer:test01

# IMPORTANT: RUN devtools::build() to build an archive, then move it to
# inst/shiny/zip and change the appropriate lines below!!!!!!!!

# get shiny plus tidyverse packages image
FROM rocker/shiny-verse:4.1.2

# install netCDF, mzR prerequisite

RUN sudo apt-get update
RUN sudo apt-get -y install netcdf-bin

# install R packages from Bioconductor

RUN R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')"

RUN R -e "BiocManager::install('mzR')"
RUN R -e "BiocManager::install('MSnbase')"

# install R packages from CRAN

RUN R -e "install.packages('enviPat', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('reshape2', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('progressr', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('magrittr', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('markdown', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('assertthat', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('glue', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('tippy', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('shinydashboard', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('shinyWidgets', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('shinyBS', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('shinyjs', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('readxl', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('Peptides', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('OrgMassSpecR', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('fs', repos='http://cran.rstudio.com/')"

# install abundancer R package from zip

COPY /zip/ /tmp/
RUN R -e "install.packages('tmp/abundancer_0.2.0.tar.gz')"

# copy the app to the image
COPY ui.R /var/abundancer/
COPY server.R /var/abundancer/
COPY about.md /var/abundancer/
COPY /www/ /var/abundancer/www/

# select port
EXPOSE 3838

# run app
CMD R -e 'shiny::runApp("/var/abundancer", port = 3838, host = "0.0.0.0")'
