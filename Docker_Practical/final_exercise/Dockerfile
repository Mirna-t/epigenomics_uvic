#!/usr/bin/env python3

FROM python:3.12-slim

LABEL maintainer="Mirna <mirna.tarzi@gmail.com>" \
      version="0.1" \
      description="Docker container for seqClass.py"


# Copy script into a directory already in PATH
COPY seqClass.py /usr/local/bin/seqClass.py

# Make it executable
RUN chmod +x /usr/local/bin/seqClass.py

# Ensure the directory is in PATH
ENV PATH="/usr/local/bin:${PATH}"
