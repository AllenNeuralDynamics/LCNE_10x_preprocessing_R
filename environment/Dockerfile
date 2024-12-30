ARG REGISTRY_HOST
FROM $REGISTRY_HOST/codeocean/py-r:python3.12.4-R4.4.2-1.2204.0-JupyterLab4.3.0-RStudiorstudio-server-2024.09.0-375-ubuntu22.04

ARG DEBIAN_FRONTEND=noninteractive

ARG GIT_ASKPASS
ARG GIT_ACCESS_TOKEN
COPY git-askpass /

ARG RCACHE_OPTIONS
RUN echo $RCACHE_OPTIONS >> $(Rscript -e 'cat(R.home())')/etc/Rprofile.site

# Set environment variables to prevent interactive prompts during installation
ENV DEBIAN_FRONTEND=noninteractive

# Update system and install required dependencies
RUN apt-get update && apt-get install -y \
    software-properties-common \
    dirmngr \
    gnupg \
    ca-certificates \
    wget \
    locales \
    --no-install-recommends \
 && locale-gen en_US.UTF-8 \
 && update-locale LANG=en_US.UTF-8

# Set the locale environment variables
ENV LANG=en_US.UTF-8 \
    LC_ALL=en_US.UTF-8

# Add the CRAN repository and GPG key for R
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 51716619E084DAB9 \
 && add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'

# Update the package index and install R
RUN apt-get update && apt-get install -y \
    r-base \
    r-base-dev \
    --no-install-recommends \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

# Optional: Install common R packages used in data science
RUN R -e "install.packages(c('tidyverse', 'data.table', 'shiny', 'rmarkdown', 'ggplot2', 'dplyr'), repos='https://cloud.r-project.org/')"

# Set the working directory
WORKDIR /workspace

# Define the command to run R interactively
CMD ["R"]
