FROM bioconductor/bioconductor_docker:RELEASE_3_16

# location for mounting
#RUN mkdir /home/rstudio/data


##----------------------------------------------------------------------------##
##  R PACKAGES
##----------------------------------------------------------------------------##

#### General CRAN Packages ####
RUN R -e 'install.packages(c("tidyverse", "remotes",  "BiocManager", "devtools", "knitr", "lintr", "styler", "here"))'

#### Project Specific CRAN Packages ####
RUN R -e 'install.packages(c("viridis", "gprofiler2", "readxl", "openxlsx"))'

#### Project Specific GitHub Packages ####
RUN R -e 'devtools::install_github("kvittingseerup/IsoformSwitchAnalyzeR")'

#### Project Specific Bioconductor Packages ####
RUN R -e 'BiocManager::install(c("ComplexHeatmap", "AnnotationDbi", "org.Mm.eg.db", "org.Hs.eg.db"))'

##----------------------------------------------------------------------------##
## SYSTEM PACKAGES
##----------------------------------------------------------------------------##

RUN apt-get update --yes && \
  apt-get upgrade --yes

RUN apt-get install --yes \
  python-dev-is-python3 \
  python3-pip \
  python3-dev \
  python3-pip \
  python3-tk

RUN apt-get install --yes \
  libcurl4-openssl-dev \
  libssl-dev \
  libcairo2-dev libxt-dev \
  libxml2-dev \
  libudunits2-dev \
  libhdf5-dev \
  libv8-dev \
  libgdal-dev \
  xorg libx11-dev libglu1-mesa-dev \
  samtools libboost-all-dev \
  libgsl-dev

##----------------------------------------------------------------------------##
## PYTHON PACKAGES
##----------------------------------------------------------------------------##


# RUN pip3 install #pythonPackageHere



