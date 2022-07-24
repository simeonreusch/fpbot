FROM python:3.10.5-buster

RUN apt-get update && apt-get install -y \
    git \
    python3-pip \
    libpq-dev \
    python3-dev \
    poetry

RUN git clone https://github.com/simeonreusch/fpbot.git
WORKDIR fpbot
EXPOSE 8000
RUN poetry install
#CMD ["mkdir", "/ZTFDATA"]
ADD .ztfquery /root
ENV ZTFDATA=/ZTFDATA/
WORKDIR /fpbot/fpbot
