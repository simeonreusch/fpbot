FROM python:3.8.5-buster

RUN apt-get update && apt-get install -y \
    git \
    python3-pip \
    libpq-dev \
    python3-dev

RUN git clone https://github.com/simeonreusch/ztffps.git
WORKDIR ztffps
EXPOSE 8000
RUN pip3 install -r requirements.txt
#CMD ["mkdir", "/ZTFDATA"]
ADD .ztfquery /root
ENV ZTFDATA=/ZTFDATA/
WORKDIR /ztffps/ztffps
