FROM ubuntu:22.04

WORKDIR /app

# Install dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    git \
    vim \
    wget

# Install Miniconda
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh
ENV PATH="/opt/conda/bin:$PATH"
RUN conda init bash

# Major conda install
COPY production-environment.yml .
RUN conda env create -f production-environment.yml

# Activate the environment
RUN echo "conda activate aceretro-env" >> ~/.bashrc
ENV PATH /opt/conda/envs/aceretro-env/bin:$PATH

# Minor installs and a huge file copy to retail model files (.pt, pytorch )
RUN pip install opennmt-py==2.3.0
# Copy just a sub-folder to optimize docker layers cache
COPY ./pathway_search_standalone/rxn_cluster_token_prompt /app/pathway_search_standalone/rxn_cluster_token_prompt
RUN pip install -e pathway_search_standalone/rxn_cluster_token_prompt

COPY . /app

# Remove the frequently changing files if they were copied
RUN rm -f /app/test_docker_build.py /app/entrypoint.py

# add frequenty changing files here (AND IN THE .dockerignore) to preserve build cache
COPY --chown=root:root test_docker_build.py /app/
COPY --chown=root:root entrypoint.py /app/

# A rigerous test of full functinoality. (OPTIONAL, Compute expensive, able to be GPU-accelerated.)
# RUN python test_docker_build.py

# default command
CMD ["/bin/bash"]
