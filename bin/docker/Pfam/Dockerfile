FROM bioconductor/bioconductor_docker:RELEASE_3_14

# location for mounting
#RUN mkdir /home/rstudio/data

##----------------------------------------------------------------------------##
## SYSTEM PACKAGES
##----------------------------------------------------------------------------##

RUN apt-get update --yes && \
  apt-get upgrade --yes

RUN apt-get install --yes \
  build-essential

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
  libgsl-dev \
  tar \
  wget \
  make 

RUN export LC_ALL=en_US.UTF-8
RUN export LANG=en_US.UTF-8
RUN locale-gen en_US.UTF-8

##----------------------------------------------------------------------------##
## Install perl
##----------------------------------------------------------------------------##

ADD http://www.cpan.org/src/5.0/perl-5.36.0.tar.gz /opt
WORKDIR /opt
RUN tar -xzf perl-5.36.0.tar.gz

WORKDIR /opt/perl-5.36.0
RUN ./Configure -des -Dprefix=/opt/perl-5.36.0/localperl

RUN make
# RUN make test
RUN make install

# Set up cpan minus
RUN cpan App::cpanminus

##----------------------------------------------------------------------------##
## Set up Pfam dependencies and data files
##----------------------------------------------------------------------------##

WORKDIR /opt
# Get PfamScan Script
RUN wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/Tools/PfamScan.tar.gz && \
  tar zxvf PfamScan.tar.gz && \
  cd PfamScan

# Get Hmmer
RUN wget http://eddylab.org/software/hmmer/hmmer.tar.gz && \
  tar zxf hmmer.tar.gz && \
  cd hmmer-3.3.2 && \
  ./configure --prefix /opt/ && \
  make && \
  make check  && \
  make install  && \
  (cd easel; make install)

# Add Hmmer to PATH
ENV PATH /opt/bin:$PATH

# Install moose dependency
RUN cpanm Moose 

# Install required bioperl pfam module and move to library
RUN cpan Bio::Pfam::Scan::PfamScan module
RUN cp -r /opt/PfamScan/Bio /usr/local/lib/x86_64-linux-gnu/perl/5.30.0/

# Install IPC module
RUN cpan IPC::Run module

# Install Pfam data files
WORKDIR /PfamScan/PfamFiles
RUN wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz && \
  wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz && \
  wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/active_site.dat.gz
  
RUN gunzip *

RUN hmmpress Pfam-A.hmm

##----------------------------------------------------------------------------##
## Run PfamScan script
##----------------------------------------------------------------------------##
WORKDIR ../../opt/PfamScan/
