FROM --platform=linux/amd64 ubuntu:20.04


RUN apt-get update && apt-get install -y ncbi-blast+ python3 python3-pip wget mafft \
&& wget https://github.com/epiliper/vapid4/releases/download/v1.0.0/vapid4_v1.0.0.tar.gz 

RUN tar -xvzf vapid4_v1.0.0.tar.gz

WORKDIR vapid4_v1.0.0

COPY requirements.txt .

RUN pip install -r requirements.txt \
&& wget https://ftp.ncbi.nlm.nih.gov/asn1-converters/by_program/table2asn/linux64.table2asn.gz \
&& gzip -d linux64.table2asn.gz \
&& mv linux64.table2asn /usr/local/bin/table2asn \
&& chmod +x /usr/local/bin/table2asn

ENV PYTHONPATH="{PYTHONPATH}:/vapid4_v1.0.0/"









