FROM rocker/tidyverse:latest

LABEL org.label-schema.license="GPL-2.0" \
      org.label-schema.vcs-url="https://github.com/covidestim/covidestim" \
      org.label-schema.vendor="Covidestim" \
      maintainer="Marcus Russi <marcus.russi@yale.edu>"

# Install node because a dependency of Covidestim has to link to it. Then,
# remove the cached repositories so they don't take up space. Without doing
# this, V8 will be unable to link to libnode.so.
RUN wget -O - https://deb.nodesource.com/setup_14.x | bash - \
        && apt-get install -y libnode-dev \
        && rm -rf /var/lib/apt/lists/*

# All future commands are run as 'rstudio' user
USER rstudio

# Copy over the GitHub repo
COPY --chown=rstudio . /tmp/covidestim-install/

# Enable O3 compilation
RUN Rscript /tmp/covidestim-install/O3-enable.R 

RUN r -e "remotes::install_deps('/tmp/covidestim-install')"

# Now install covidestim
# RUN r -e "devtools::install('/tmp/covidestim-install')" \
#   && rm -rf /tmp/covidestim-install

RUN R CMD INSTALL --preclean "/tmp/covidestim-install" \
  && rm -rf /tmp/covidestim-install

CMD ["R"]
