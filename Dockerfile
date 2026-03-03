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
# RUN apt-get update && apt-get install -y --no-install-recommends \
#     build-essential \
#     ca-certificates \
#     curl \
#     gcc \
#     g++ \
#     git \
#     libxext6 \
#     make \
#     wget \
#     && update-ca-certificates \
#     && rm -rf /var/lib/apt/lists/*
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential git cmake wget libboost-all-dev libeigen3-dev \
    libgoogle-glog-dev libprotobuf-dev protobuf-compiler libhdf5-dev \
    libatlas-base-dev python3-dev librdkit-dev python3-numpy python3-pip python3-pytest \
    libjsoncpp-dev \
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

USER root
# RUN git clone https://github.com/dkoes/openbabel.git && \
#     cd openbabel && \
#     mkdir build && \
#     cd build && \
#     cmake -DWITH_MAEPARSER=OFF -DWITH_COORDGEN=OFF -DPYTHON_BINDINGS=ON -DRUN_SWIG=ON .. && \
#     make && \
#     make install
# RUN git clone https://github.com/gnina/libmolgrid.git && \
#     cd libmolgrid && \
#     mkdir build && \
#     cd build && \
#     cmake -DOPENBABEL3_INCLUDE_DIR=$MAMBA_ROOT_PREFIX/envs/env/include .. && \
#     make && \
#     make install
RUN git clone https://github.com/gnina/gnina.git && \
    cd gnina && \
    mkdir build && \
    cd build && \
    cmake ..  \
        -DCUDAToolkit_ROOT=$MAMBA_ROOT_PREFIX/envs/env \
        -DLIBMOLGRID_INCLUDE_DIR=$MAMBA_ROOT_PREFIX/envs/env/include \
        -DLIBMOLGRID_LIBRARY=$MAMBA_ROOT_PREFIX/envs/env/lib \
        -DZLIB_ROOT=$MAMBA_ROOT_PREFIX/envs/env \
        -DZLIB_LIBRARY=$MAMBA_ROOT_PREFIX/envs/env/lib/libz.so \
        -DZLIB_INCLUDE_DIR=$MAMBA_ROOT_PREFIX/envs/env/include && \
    make && \
    make install
# -DUSE_SYSTEM_NVTX=1 may be needed with pytorch 2.7.0 and CUDA 12.9
# RUN git clone https://github.com/gnina/gnina.git && \
#     cd gnina && \
#     mkdir build && \
#     cd build && \
#     cmake .. \
#         -DCMAKE_CUDA_ARCHITECTURES=all \
#         -DUSE_SYSTEM_NVTX=1 \
#         -DCUDAToolkit_ROOT=$MAMBA_ROOT_PREFIX/envs/env \
#         -DCMAKE_PREFIX_PATH="$(python -c 'import torch; print(torch.utils.cmake_prefix_path)');$MAMBA_ROOT_PREFIX/envs/env" \
#         -DNumPy_DIR=$MAMBA_ROOT_PREFIX/envs/env/lib/python3.??/site-packages/numpy/_core/include \
#         -DCMAKE_INSTALL_PREFIX=$MAMBA_ROOT_PREFIX/envs/env \
#         -DOPENBABEL3_INCLUDE_DIR=$MAMBA_ROOT_PREFIX/envs/env/include \
#         -DLIBMOLGRID_INCLUDE_DIR=$MAMBA_ROOT_PREFIX/envs/env/include \
#         -DLIBMOLGRID_LIBRARY=$MAMBA_ROOT_PREFIX/envs/env/lib \
#         -DRDKIT_INCLUDE_DIR=$MAMBA_ROOT_PREFIX/envs/env/include \
#         -DRDKIT_LIBRARIES=$MAMBA_ROOT_PREFIX/envs/env/lib \
#         -DZLIB_ROOT=$MAMBA_ROOT_PREFIX/envs/env \
#         -DZLIB_LIBRARY=$MAMBA_ROOT_PREFIX/envs/env/lib/libz.so \
#         -DZLIB_INCLUDE_DIR=$MAMBA_ROOT_PREFIX/envs/env/include && \
#     make && \
#     make install
USER 1000

# Smoke test: confirm key tools are on PATH and importable
RUN gnina --help | head -1 && \
    fpocket -h 2>&1 | head -1 && \
    python -c "from rdkit import Chem; print('RDKit', Chem.rdBase.rdkitVersion)"
