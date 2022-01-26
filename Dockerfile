FROM continuumio/miniconda3

#Activate Shell
SHELL ["/bin/bash", "-c"]
ENV PATH="/opt/conda/bin/:$PATH"

RUN apt-get update
RUN apt-get -y install gcc
# RUN apt-get -y install python3.6

RUN conda update -n base -c defaults conda
RUN conda create -n htpmd python=3.6

WORKDIR /src

ENV PATH="/opt/conda/envs/htpmd/bin:$PATH"

COPY /src .

RUN source /opt/conda/bin/activate htpmd && \
    conda install pytorch==1.8.1 cpuonly -c pytorch && \
    conda install pyg -c pyg -c conda-forge && \
    conda install -c conda-forge mdtraj && \
    pip install --upgrade pip && \
    pip install -e .[tests] && \
    pip install -e .
