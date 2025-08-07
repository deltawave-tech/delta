FROM debian:stable

LABEL maintainer="PharmAI GmbH <contact@pharm.ai>" \
    org.label-schema.name="PLIP: The Protein-Ligand Interaction Profiler" \
    org.label-schema.description="https://www.doi.org/10.1093/nar/gkv315"

ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies including PLIP requirements
RUN apt-get update && apt-get install -y \
    pymol \
    python3-distutils \
    python3-lxml \
    python3-pymol \
    python3-numpy \
    python3-pip \
    python3-venv \
    python3-dev \
    python3-openbabel \
    openbabel \
    build-essential; \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Setup virtual environment with access to system packages
ENV VIRTUAL_ENV=/opt/venv
RUN python3 -m venv $VIRTUAL_ENV --system-site-packages
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

# Install API dependencies
COPY requirements.txt /tmp/
RUN pip install --no-cache-dir -r /tmp/requirements.txt

# copy PLIP source code
WORKDIR /src
COPY . .
RUN chmod +x plip/plipcmd.py
ENV PYTHONPATH=/src

# Create storage directory for results
RUN mkdir -p /storage && chmod 777 /storage

# Switch entry point to API
EXPOSE 8000
CMD ["python", "-m", "uvicorn", "plip.plip_api:app", "--host", "0.0.0.0", "--port", "8000", "--reload", "--reload-include", "*.py"]
