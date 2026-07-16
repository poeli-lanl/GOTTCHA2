# Stage 1: Builder
FROM mambaorg/micromamba:latest AS builder

USER root

# Install build dependencies
RUN apt-get update && \
    apt-get install -y ca-certificates && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

USER $MAMBA_USER

# Copy environment file and create conda environment
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml


RUN micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

# Copy project files and install
WORKDIR /app
COPY --chown=$MAMBA_USER:$MAMBA_USER . /app

RUN micromamba run -n base pip install --no-cache-dir .

# Stage 2: Runtime
FROM mambaorg/micromamba:latest

USER root

# Install minimal runtime dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        ca-certificates && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

USER $MAMBA_USER

# Copy conda environment from builder
COPY --from=builder /opt/conda /opt/conda

# Set working directory
WORKDIR /data

# Set environment
ENV PATH="/opt/conda/bin:$PATH" \
    PYTHONUNBUFFERED=1

# Verify installation
RUN gottcha2 version 2>/dev/null || echo "GOTTCHA2 installed"

# Default command
CMD ["gottcha2", "--help"]
