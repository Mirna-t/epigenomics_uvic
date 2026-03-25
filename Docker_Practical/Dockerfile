#!/usr/bin/env Rscript

# Start from debian linux image (DockerHub)
FROM debian:stable

# Add custom label
LABEL maintainer="Mirna <mirna.tarzi@gmail.com>" \
      version="0.1" \
      description="Docker container for Docker practical"

# Install R (after apt-get update)
RUN apt-get update && apt-get install -y r-base

# Install R package 'optparse'
RUN R -e 'install.packages("optparse", repos = "http://cloud.r-project.org/")'

# Make the folder '/scripts' in the container
RUN mkdir /scripts

# Copy 'scripts/script.R' to the folder '/scripts' in the container
ADD scripts/script.R /scripts

# Give execution permissions to script.R
RUN chmod +x /scripts/script.R

# Add /scripts folder to the PATH environment variable
ENV PATH="$PATH:/scripts"           
