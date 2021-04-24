FROM python:3.8-slim
WORKDIR /home/biolib
COPY requirements.txt .
RUN pip install -r requirements.txt
RUN apt-get install dssp
RUN apt-get install build-essentials wget libjson0 libjson0-dev libxml2 libxml2-dev libxml2-utils autoconf automake libtool
COPY predict.py .
COPY data/test.zip data/
ENTRYPOINT ["python3", "predict.py"]
