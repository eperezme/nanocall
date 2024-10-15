# Use Ubuntu as the base image
FROM nvidia/cuda:11.7.1-runtime-ubuntu20.04

# Set non-interactive mode for apt
ENV DEBIAN_FRONTEND=noninteractive

# Update and install dependencies
RUN apt-get update && \
    apt-get install -y \
    wget \
    git \
    curl \
    ca-certificates \
    software-properties-common \
    build-essential \
    cmake \
    python3 \
    python3-pip \
    libncurses5-dev \
    libncursesw5-dev \
    zlib1g-dev && \
    rm -rf /var/lib/apt/lists/*

# Install ONT Dorado
WORKDIR /opt
RUN wget https://cdn.oxfordnanoportal.com/software/analysis/dorado-0.8.1-linux-x64.tar.gz && \
    tar -xzf dorado-0.8.1-linux-x64.tar.gz && \
    rm dorado-0.8.1-linux-x64.tar.gz && \
    mv dorado-0.8.1-linux-x64 dorado

# Set PATH for Dorado
ENV PATH="/opt/dorado/bin:${PATH}"
RUN mkdir -p /opt/dorado/models
ENV DORADO_MODELS_DIRECTORY="/opt/dorado/models"

# Define a build argument to control model download
ARG DOWNLOAD_MODELS=false

# Conditionally download models based on the build argument
RUN if [ "$DOWNLOAD_MODELS" = "true" ]; then \
        dorado download --model all; \
    else \
        echo "Skipping model download"; \
    fi

# Set up entrypoint for running ONT Dorado
WORKDIR /data

CMD ["/bin/bash"]
