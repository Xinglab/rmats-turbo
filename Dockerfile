#rmats dockerfile without pairadise
# usage: python /usr/bin/rmats.py -h
FROM conda/miniconda3
LABEL maintainer="sridhar sridhar@wustl.edu"
LABEL docker_image rmats_turbo4.1.0

# Fundamentals

RUN apt-get update -y && apt-get install -y --no-install-recommends \
    build-essential \
    bzip2 \
    curl \
    cmake \
    g++ \
    git \
    gcc \
    gfortran \ 
    less \
    libcurl4-openssl-dev \
    libpng-dev \
    libopencv-dev \
    libgsl-dev \
    libblas-dev \
    liblapack-dev \
    vim-tiny \
    libnss-sss \
    libssl-dev \
    libxml2-dev \
    make \
    pkg-config \
    rsync \
    unzip \
    wget \
    zip \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    hdf5-tools \
    libhdf5-dev \
    hdf5-helpers \
    ncurses-dev
    

RUN pip install --upgrade pip && \
    pip install setuptools && \
    pip install numpy && \
    pip install matplotlib && \
    pip install pandas && \
    pip install scipy && \
    pip install pysam && \
    pip install biopython && \
    pip install seaborn && \
    pip install scikit-learn && \
    pip install tqdm && \
    pip install --upgrade cython
    
#Create Working Directory
WORKDIR /docker_main
RUN wget https://github.com/Xinglab/rmats-turbo/releases/download/v4.1.0/rmats_turbo_v4_1_0.tar.gz && \
    tar -xzf rmats_turbo_v4_1_0.tar.gz && \
    cd rmats-turbo && \
    cp -r * /usr/bin && \
    chmod +x /usr/bin/build_rmats && \
    /usr/bin/build_rmats --no-paired-model && \
    chmod +x /usr/bin/rmats.py



#install samtools
WORKDIR /docker_main
RUN wget https://github.com/samtools/samtools/releases/download/1.4/samtools-1.4.tar.bz2 && \
    tar -jxf samtools-1.4.tar.bz2 && \
    cd samtools-1.4 && \
    make && \
    make install && \
    cp samtools /usr/bin/

# Clean up
RUN cd /docker_main / && \
   rm -rf samtools-1.4 rmats_turbo_v4_1_0.tar.gz && \
   apt-get autoremove -y && \
   apt-get autoclean -y  && \
   apt-get clean && \
   apt-get clean all && rm -rf /var/lib/apt/lists/*

# Set default working path
WORKDIR /docker_main
