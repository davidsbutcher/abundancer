# DOCKERFILE for ABUNDANCER
# github.com/davidsbutcher/abundancer
# docker build --tag davidsbutcher/abundancer:test01 --no-cache .

# get shiny plus tidyverse packages image
FROM rocker/shiny-verse:4.0.3

# install R packages from Bioconductor

RUN R -e "if (!requireNamespace('BiocManager', quietly = TRUE)) install.packages('BiocManager')"

RUN R -e "BiocManager::install('mzR')"
RUN R -e "BiocManager::install('MSnbase')"

# install R packages from github

RUN R -e "install.packages('devtools', repos='http://cran.rstudio.com/')"

RUN R -e "devtools::install_github('davidsbutcher/abundancer')"

# install R packages from CRAN 

RUN R -e "install.packages('magrittr', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('assertthat', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('glue', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('tippy', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('shinydashboard', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('shinyWidgets', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('shinyjs', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('readxl', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('Peptides', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('OrgMassSpecR', repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('fs', repos='http://cran.rstudio.com/')"

# copy the app to the image
COPY ui.R /var/abundancer/
COPY server.R /var/abundancer/
COPY about.md /var/abundancer/
COPY /www/ /var/abundancer/www/

# select port
EXPOSE 3838

# run app
CMD R -e 'shiny::runApp("/var/abundancer", port = 3838, host = "0.0.0.0")'