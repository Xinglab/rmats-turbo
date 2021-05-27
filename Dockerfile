FROM debian:buster

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
       ca-certificates \
       cmake \
       curl \
       cython \
       g++ \
       gfortran \
       git \
       libblas-dev \
       libgsl-dev \
       liblapack-dev \
       make \
       python-dev \
       python-numpy \
       r-base \
       r-cran-nloptr \
       zlib1g-dev \
    && rm -rf /var/lib/apt/lists/* \
    # Use a build dir to be removed after artifacts are extracted
    && mkdir /rmats_build \
    && cd /rmats_build \
    && git clone https://github.com/Xinglab/rmats-turbo.git \
    && cd rmats-turbo \
    # The build will source setup_environment.sh which will source ~/.bashrc.
    # Skip that by truncating setup_environment.sh
    && echo '' > setup_environment.sh \
    && ./build_rmats \
    # Copy the build results
    && mkdir /rmats \
    && cd /rmats \
    && cp /rmats_build/rmats-turbo/rmats.py ./ \
    && cp /rmats_build/rmats-turbo/cp_with_prefix.py ./ \
    && cp /rmats_build/rmats-turbo/*.so ./ \
    && mkdir rMATS_C \
    && cp /rmats_build/rmats-turbo/rMATS_C/rMATSexe ./rMATS_C \
    && mkdir rMATS_P \
    && cp /rmats_build/rmats-turbo/rMATS_P/*.py ./rMATS_P \
    && mkdir rMATS_R \
    && cp /rmats_build/rmats-turbo/rMATS_R/*.R ./rMATS_R \
    # Remove build dir
    && rm -rf /rmats_build \
    # Build STAR
    && mkdir /star_build \
    && cd /star_build \
    && curl -L -O https://github.com/alexdobin/STAR/archive/refs/tags/2.7.9a.tar.gz \
    && tar -xvf 2.7.9a.tar.gz \
    && cd STAR-2.7.9a/source \
    && make STAR \
    && cp STAR /usr/local/bin

# Set defaults for running the image
WORKDIR /rmats
ENTRYPOINT ["python", "/rmats/rmats.py"]
CMD ["--help"]
