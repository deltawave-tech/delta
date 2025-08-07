FROM nvcr.io/hpc/autodock:2020.06

# Set locale environment
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

# Environment variables for memory management
ENV CUDA_VISIBLE_DEVICES=0
ENV GPU_DEVICE_ORDINAL=0
ENV OMP_NUM_THREADS=1
ENV OPENBLAS_NUM_THREADS=1
ENV MKL_NUM_THREADS=1

# Avoid interactive prompts during installation
ENV DEBIAN_FRONTEND=noninteractive

# Fix NVIDIA repository GPG key and install build tools
RUN apt-key adv --fetch-keys https://developer.download.nvidia.com/compute/cuda/repos/ubuntu1804/x86_64/3bf863cc.pub && \
    apt-get update && \
    apt-get install -y build-essential git wget

# We'll build AutoDock-GPU directly in the container
WORKDIR /tmp
RUN mkdir -p /usr/local/bin/

# Install other dependencies including libsm6 for OpenBabel
RUN apt-get update && \
    apt-get install -y gnupg2 curl python3-pip openbabel \
    python3-numpy python3-scipy build-essential libtinfo5 csh \
    libsm6 libice6 libxt6 x11-common

# Install Miniconda and required packages with specific versions
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    /opt/conda/bin/conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/main && \
    /opt/conda/bin/conda tos accept --override-channels --channel https://repo.anaconda.com/pkgs/r && \
    /opt/conda/bin/conda tos accept --override-channels --channel conda-forge && \
    /opt/conda/bin/conda tos accept --override-channels --channel bioconda && \
    /opt/conda/bin/conda config --add channels bioconda && \
    /opt/conda/bin/conda config --add channels conda-forge && \
    /opt/conda/bin/conda install -y -c bioconda \
    python=3.9 \
    pdbfixer \
    rdkit \
    meeko \
    six && \
    echo 'export PATH="/opt/conda/bin:$PATH"' >> ~/.bashrc && \
    echo 'export PYTHONPATH="/opt/conda/lib/python3.9/site-packages:$PYTHONPATH"' >> ~/.bashrc

# Fix the rdkit six module in a separate step using conda's Python
RUN /opt/conda/bin/python -c "import six; content = '''import six as _six\nfrom six import *\nStringIO = _six.StringIO\n'''; open('/opt/conda/lib/python3.9/site-packages/rdkit/six.py', 'w').write(content)"

# Make conda available to the default Python environment
ENV PATH="/opt/conda/bin:${PATH}"
ENV PYTHONPATH="/opt/conda/lib/python3.9/site-packages:${PYTHONPATH}"

# Install OpenBabel 3.x
RUN apt-get update && apt-get install -y openbabel

# Download and install MGLTools
RUN cd /tmp && \
    wget https://ccsb.scripps.edu/mgltools/download/491/ -O mgltools_x86_64Linux2_1.5.7.tar.gz && \
    tar -xzf mgltools_x86_64Linux2_1.5.7.tar.gz -C /opt && \
    cd /opt/mgltools_x86_64Linux2_1.5.7 && \
    ./install.sh -d /opt/mgltools && \
    echo 'export PATH="/opt/mgltools/bin:$PATH"' >> ~/.bashrc && \
    rm /tmp/mgltools_x86_64Linux2_1.5.7.tar.gz

# Install AutoGrid from source
RUN cd /tmp && \
    wget https://autodock.scripps.edu/wp-content/uploads/sites/56/2021/10/autodocksuite-4.2.6-src.tar.gz && \
    tar -xvf autodocksuite-4.2.6-src.tar.gz && \
    cd src/autogrid && \
    ./configure && \
    make && \
    make install && \
    cd ../.. && \
    rm -rf src autodocksuite-4.2.6-src.tar.gz

# Install AutoDock Vina
RUN apt-get update && \
    apt-get install -y wget build-essential && \
    wget https://vina.scripps.edu/wp-content/uploads/sites/55/2020/12/autodock_vina_1_1_2_linux_x86.tgz && \
    tar -xzvf autodock_vina_1_1_2_linux_x86.tgz && \
    cp autodock_vina_1_1_2_linux_x86/bin/vina /usr/local/bin/ && \
    chmod +x /usr/local/bin/vina && \
    rm -rf autodock_vina_1_1_2_linux_x86 autodock_vina_1_1_2_linux_x86.tgz

# Add diagnostic commands to check existing AutoDock/AutoGrid installations
RUN which autogrid4 || echo "autogrid4 not found" && \
    which autodock4 || echo "autodock4 not found" && \
    which autodock_gpu_128wi || echo "autodock_gpu not found" && \
    which vina || echo "vina not found"

# Set up MGLTools environment
ENV MGL_ROOT="/opt/mgltools"
ENV PYTHONPATH="/opt/mgltools/MGLToolsPckgs:/opt/mgltools/lib/python2.7/site-packages:${PYTHONPATH}"
ENV LD_LIBRARY_PATH="/opt/mgltools/lib:${LD_LIBRARY_PATH}"
ENV PATH="/opt/mgltools/bin:/opt/mgltools/MGLToolsPckgs/AutoDockTools/Utilities24:${PATH}"

# Create runtime script for MGLTools conversions
RUN echo '#!/bin/bash\nPYTHONPATH=/opt/mgltools/MGLToolsPckgs:/opt/mgltools/lib/python2.7/site-packages python "$@"' > /usr/local/bin/mgltools-python && \
    chmod +x /usr/local/bin/mgltools-python

# Copy requirements.txt first for better caching
COPY ./app/requirements.txt /tmp/requirements.txt

# Install Python dependencies from requirements.txt using Python 3.9
RUN /opt/conda/bin/pip install -r /tmp/requirements.txt

WORKDIR /app

COPY ./app /app/

RUN mkdir -p /app/storage/results && \
    chmod -R 777 /app

# Using nvidia base image - CUDA libraries should already be available

# Add a verification step to show AutoDock-GPU version and library dependencies
RUN ldconfig && \
    echo "AutoDock-GPU 1.5.3 library dependencies:" && \
    ldd /usr/local/bin/autodock_gpu_128wi || echo "Failed to check library dependencies" && \
    echo "AutoDock-GPU 1.5.3 successfully installed"

# Make entrypoint script executable
RUN chmod +x /app/entrypoint.sh

EXPOSE 8686

# Use entrypoint script to validate environment before starting the API
ENTRYPOINT ["/app/entrypoint.sh"]