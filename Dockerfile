FROM continuumio/miniconda3

WORKDIR /home/htp_md
COPY ./ ./

#Activate Shell
SHELL ["/bin/bash", "-c"]
ENV PATH="/opt/conda/bin/:$PATH"
RUN conda update -n base -c defaults conda
RUN conda env create -f env.yml

WORKDIR /home/htp_md
ENV PATH="/opt/conda/envs/md_worker/bin:$PATH"

COPY . /home/htp_md

RUN source /opt/conda/bin/activate htpmd && \
    pip install -e .[tests]

WORKDIR /home/htp_md
RUN python setup.py install

CMD ["tail", "-f", "/dev/null"]
