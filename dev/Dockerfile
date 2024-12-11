# Author: Etienne CAMENEN
# Date: 2024
# Contact: etienne.camenen@gmail.com

FROM rocker/rstudio:4.0

MAINTAINER Etienne CAMENEN (etienne.camenen@gmail.com)

ENV _R_CHECK_FORCE_SUGGESTS_ FALSE
ENV PKGS cmake git libcurl4-openssl-dev libglu1-mesa liblapack-dev libssl-dev libxml2-dev qpdf texlive-fonts-recommended texlive-latex-extra texlive-latex-recommended
ARG TOOL_NAME
ARG TOOL_VERSION

RUN apt-get update --allow-releaseinfo-change -qq && \
    apt-get install -y ${PKGS}
RUN apt-get install -y --no-install-recommends libglpk-dev libxt6
ENV RPKGS BiocManager config covr cowplot DataExplorer devtools janitor knitr lintr magrittr markdown openxlsx plotly RColorBrewer reshape2 rlist rmarkdown rstatix tidyverse
RUN Rscript -e "install.packages(commandArgs(TRUE))" ${RPKGS}
RUN Rscript -e "BiocManager::install('BiocCheck')"
RUN R -e "devtools::install_github('ecamenen/"${TOOL_NAME}"', ref = '"${TOOL_VERSION}"')"
RUN apt-get purge -y git g++ && \
	apt-get autoremove --purge -y && \
	apt-get clean && \
	rm -rf /var/lib/{cache,log}/ /tmp/* /var/tmp/*
COPY .. /home/rstudio
