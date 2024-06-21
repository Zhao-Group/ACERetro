FROM ubuntu:22.04

WORKDIR /app

# Install dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    git \
    wget

# Install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh

ENV PATH="/opt/conda/bin:$PATH"

# Major conda install
COPY kas-environment.yml .
RUN conda env create -f kas-environment.yml

# Activate the environment
RUN echo "conda activate aceretro-from-env-v5" >> ~/.bashrc
ENV PATH /opt/conda/envs/aceretro-from-env-v5/bin:$PATH

RUN pip install opennmt-py==2.3.0

COPY . /app
RUN pip install -e pathway_search_standalone/rxn_cluster_token_prompt

# Compute expensive, able to be GPU-accelerated. It's a rigerous test of full functinoality.
RUN python test_docker_build.py
# Define default command
CMD ["python", "test_docker_build.py"]