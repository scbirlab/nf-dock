FROM mambaorg/micromamba:1.5.6

LABEL org.opencontainers.image.title="nf-dock"
LABEL org.opencontainers.image.description="GNINA docking pipeline"
LABEL org.opencontainers.image.authors="SCBIR Lab, Francis Crick Institute"

ENV HOME=/home/micromamba
ENV MAMBA_ROOT_PREFIX=/opt/conda
WORKDIR $HOME

USER root
RUN echo "user:x:1001:1001::/home/user:/bin/bash" >> /etc/passwd && \
    mkdir -p /home/user && chown -R 1001:1001 /home/user
RUN mkdir -p $HOME/.conda && chown -R 1000:1000 $HOME
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    ca-certificates \
    curl \
    gcc \
    g++ \
    && update-ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# NVIDIA runtime: tell the container to use host GPU driver
# These are picked up by nvidia-container-toolkit / Singularity --nv
ENV NVIDIA_VISIBLE_DEVICES=all
ENV NVIDIA_DRIVER_CAPABILITIES=compute,utility
USER 1000

COPY environment.yml /tmp/environment.yml
RUN micromamba create -n env -f /tmp/environment.yml && \
    micromamba clean --all --yes

ENV PATH=$MAMBA_ROOT_PREFIX/envs/env/bin:$PATH

# Smoke test: confirm key tools are on PATH and importable
RUN fpocket -h 2>&1 | head -1 && \
    python -c "from rdkit import Chem; print('RDKit', Chem.rdBase.rdkitVersion)"
