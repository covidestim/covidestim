 ## Emacs, make this -*- mode: sh; -*-

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

ENV R_BASE_VERSION 4.0.2

## Now install covidestim
RUN installGithub.r covidestim/covidestim \
	&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds

CMD ["R"]
