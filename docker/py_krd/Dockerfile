# golob/gappa:0.1

FROM --platform=amd64 ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y \
git \
python3-pip \
cmake \
&& apt-get clean \
&& apt-get purge \
&& rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN mkdir /src/

WORKDIR /src/
RUN git clone --recursive https://github.com/lczech/gappa.git && cd gappa/ && make && cp bin/* /usr/local/bin/ && cd /src/ && rm -r /src/gappa/
WORKDIR /root/

ADD interval_krd.py /usr/local/bin/
RUN chmod +x /usr/local/bin/interval_krd.py