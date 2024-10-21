# using barebones conda base
# pipeline needs to be installed with python 3
# make sure docker base is using python 3
FROM condaforge/mambaforge as ubuntu
ENV HOME /root

# Update Ubuntu
RUN apt-get update && apt update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends tzdata && apt-get -y --fix-missing install binutils build-essential wget cmake make gcc g++ bzip2 tar zlib1g-dev libjpeg-dev gfortran bison flex file fish libhe5-hdfeos-dev liblapack-dev libblas-dev libboost-dev libarmadillo-dev libmbedtls-dev

FROM ubuntu as mamba

# Prevent mamba from using the default conda channels run by Anaconda/NumFOCUS
RUN echo " \
allowlist_channels: \
  - conda-forge \
channels: \
  - conda-forge \
channel_priority: strict \
default_channels: [] \
" > /opt/conda/.condarc
# Update the base environment
RUN mamba update --all --yes

FROM mamba as python

RUN mamba install -y -c conda-forge "python=3.11"

FROM python as python_environment

COPY ECOSTRESS.yml /root/ECOSTRESS.yml
RUN mamba env update -n base -f /root/ECOSTRESS.yml

FROM python_environment as julia

# Set environment variables for Julia installation
ENV JULIA_VERSION=1.10

# Download and install Julia
RUN wget https://julialang-s3.julialang.org/bin/linux/x64/${JULIA_VERSION}/julia-${JULIA_VERSION}-latest-linux-x86_64.tar.gz && \
    tar xzf julia-${JULIA_VERSION}-latest-linux-x86_64.tar.gz -C /usr --strip-components 1 && \
    rm -rf julia-${JULIA_VERSION}-latest-linux-x86_64.tar.gz

FROM julia as julia_build

RUN mkdir /julia
COPY ./VNP43NRT_jl /julia/VNP43NRT_jl
COPY ./STARS_jl /julia/STARS_jl

RUN mkdir /julia/env
COPY ./julia_env/Project.toml /julia/env/Project.toml
COPY ./julia_env/Manifest.toml /julia/env/Manifest.toml

ENV JULIA_LOAD_PATH=/julia/env:

RUN julia -e "using Pkg; Pkg.activate(); Pkg.instantiate();"

FROM julia_build as installation

# creating directory inside container to store PGE code
RUN mkdir /app
RUN mkdir /pge
# copying current snapshot of PGE code repository into container
COPY . /app/
COPY ./PGE/*.sh /pge/
RUN chmod +x /pge/*.sh
# running install script to create the "ECOSTRESS" conda environment
WORKDIR /app

FROM installation as build
ENV CONDA_PREFIX=/opt/conda
RUN python setup.py install; rm -rvf build; rm -rvf dist; rm -rvf *.egg-info; rm -rvf CMakeFiles
