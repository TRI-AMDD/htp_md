FROM python:3.7

# RUN apt-get update
# RUN apt-get install -y git
# RUN apt-get install -y build-essential

COPY . /home/htp_md

WORKDIR /home/htp_md

RUN pip install --upgrade pip
RUN pip install -e .[tests]

CMD ["tail", "-f", "/dev/null"]
