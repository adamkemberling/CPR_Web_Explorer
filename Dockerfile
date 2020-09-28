### Dockerfile

# Docker Base Image for R-Shiny Environment
FROM rocker/shiny-verse:3.6.1

# Installing Spatial Library Dependencies
RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    libgdal-dev \
    libgeos-dev \
    libproj-dev \
    libudunits2-dev 

# Installing R Package Dependencies
RUN install2.r --error \
    tidyverse \
    ggpmisc \
    here \
    lubridate \
    mgcv \
    patchwork \
    plotly \
    rgdal \
    rgeos \
    sf \
    devtools

# Attempt at Installing gmRi Package
RUN R -e "devtools::install_github('gulfofmaine/gmri', upgrade = 'never')"

# Shiny Server Customizations.
COPY shiny-customized.config /etc/shiny-server/shiny-server.conf
COPY shiny_logs /var/log/shiny-server

# Copy the current folder into the path of the app
COPY . /srv/shiny-server

# Set working directory to the app
WORKDIR /srv/shiny-server

EXPOSE 3838

# Set the unix commands to run the app
CMD ["R", "-e","shiny::runApp('app.R', launch.browser = FALSE, port = 3838, host = '0.0.0.0')"]