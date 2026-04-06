ARG REGISTRY_HOST
FROM $REGISTRY_HOST/codeocean/py-r:python3.12.4-R4.4.2-1.2204.0-JupyterLab4.3.0-RStudiorstudio-server-2024.09.0-375-ubuntu22.04

ARG DEBIAN_FRONTEND=noninteractive


ARG GIT_ASKPASS
ARG GIT_ACCESS_TOKEN
COPY git-ask-pass /


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
 && update-locale LANG=en_US.UTF-8 \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

# Set the locale environment variables
ENV LANG=en_US.UTF-8 \
    LC_ALL=en_US.UTF-8



RUN R -q -e "install.packages(c('tidyverse', 'data.table', 'shiny', 'rmarkdown', 'ggplot2', 'dplyr'), repos='https://cloud.r-project.org/')"
RUN R -q -e "install.packages('remotes', repos='https://cloud.r-project.org/')" \
 && R -q -e "remotes::install_version('plotly',      version = '4.12.0', repos = 'https://cran.r-project.org')" \
 && R -q -e "remotes::install_version('harmony',      version = '1.2.3', repos = 'https://cran.r-project.org')" \
 && R -q -e "remotes::install_version('Seurat',       version = '5.4.0', repos = 'https://cran.r-project.org')" \
 && R -q -e "remotes::install_version('SeuratObject', version = '5.3.0', repos = 'https://cran.r-project.org')"




RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org/')"
RUN R -e "remotes::install_version('SingleCellExperiment', version = '1.28.1', repos = BiocManager::repositories())"
RUN R -e "remotes::install_version('zellkonverter', version = '1.16.0', repos = BiocManager::repositories())"

# Set the working directory
WORKDIR /workspace

# Define the command to run R interactively
CMD [\"R\"]