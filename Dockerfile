# Building the image with docker slim. See my corresponding blog post: https://kastanday.com/docker-slim
# First build the base image sudo docker build -t kastanday/aceretro:prod . 
# then slim it down (make sure to add the 4 MINIO env vars):
# sudo docker run -it --rm -v /var/run/docker.sock:/var/run/docker.sock \
#   dslim/slim build \
#   --target kastanday/aceretro:prod \
#   --tag kastanday/aceretro:slimer \
#   --http-probe=false \
#   --include-path=/opt/conda/envs/aceretro-env/lib/python3.6 \
#   --include-path=/app/pathway_search_standalone/rxn_cluster_token_prompt/rxn_cluster_token_prompt/ \
#   --include-path=/app/pathway_search_standalone/askcos-core/askcos/__init__.py \
#   --include-path=/app/pathway_search_standalone/askcos-core/askcos/global_config.py \
#   --include-path=/app/pathway_search_standalone/askcos-core/askcos/__pycache__ \
#   --include-path=/app/pathway_search_standalone/askcos-core/askcos/application \
#   --include-path=/app/pathway_search_standalone/askcos-core/askcos/global_config.py \
#   --include-path=/app/pathway_search_standalone/askcos-core/askcos/interfaces \
#   --include-path=/app/pathway_search_standalone/askcos-core/askcos/prioritization \
#   --include-path=/app/pathway_search_standalone/askcos-core/askcos/retrosynthetic \
#   --include-path=/app/pathway_search_standalone/askcos-core/askcos/synthetic \
#   --include-path=/app/pathway_search_standalone/askcos-core/askcos/utilities \
#   --include-path=/app/pathway_search_standalone/askcos-core/askcos/data/bkms-data/models/bkms/1 \
#   --include-path=/app/pathway_search_standalone/askcos-core/askcos/data/bkms-data/templates/bkms-and-reaxys-templates.json.gz \
#   --include-path=/app/pathway_search_standalone/askcos-core/askcos/data/models/fast_filter/1 \
#   --include-path=/app/pathway_search_standalone/scripts \
#   --include-path=/app/enzy_template_id_search.py \
#   --include-path=/app/draw_chemical_svg.py \
#   --include-path=/app/bkms-retro.templates.bkms.json.gz \
#   --include-path=/app/bkms-templates.df.json.gz \
#   --include-path=/app/process_reaction_database/saved_model/ecfp4_4096_3_layer_epoch10.pt \
#   --include-path=/app/pathway_search_standalone/rxn_cluster_token_prompt/models/class_token_retro_std_pistachio_201002_12tokens_with-reagents_step_260000.pt \
#   --include-path=/app/pathway_search_standalone/rxn_cluster_token_prompt/models/data_balancing_forward_std_pistachio_201002_with-reagents_augm-rot-rotmolorder_step_255000.pt \
#   --include-path=/app/pathway_search_standalone/rxn_cluster_token_prompt/models/data_balancing_classification_std_pistachio_201002_seq2seq_alone_step_500000.pt \
#   --include-bin=/usr/lib/x86_64-linux-gnu/libXrender.so.1 \
#   --include-path=/usr/lib/x86_64-linux-gnu \
#   --include-shell \
#   --include-oslibs-net \
#   --env MINIO_URL=XXX \
#   --env MINIO_ACCESS_KEY=XXX \
#   --env MINIO_SECRET_ACCESS_KEY=XXX \
#   --env MINIO_SECURE=False \
#   --exec "python entrypoint.py --job_id 21"

FROM ubuntu:22.04

WORKDIR /app

# Install dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    git \
    vim \
    wget \
    libxrender1

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
ENV PATH="/opt/conda/envs/aceretro-env/bin:$PATH"

# Installs that can't be in production-environment.yaml because Conda refuses to install w/ dependency conflicts.
RUN pip install opennmt-py==2.3.0
# Copy just a sub-folder to optimize docker layers cache
COPY ./pathway_search_standalone/rxn_cluster_token_prompt /app/pathway_search_standalone/rxn_cluster_token_prompt
RUN pip install -e pathway_search_standalone/rxn_cluster_token_prompt
RUN pip install rxnmapper==0.4.0

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
