FROM python:3.6

COPY . /home/htp_md

WORKDIR /home/htp_md

RUN pip install --upgrade pip
RUN pip install -e .[tests]
RUN pip install -e .