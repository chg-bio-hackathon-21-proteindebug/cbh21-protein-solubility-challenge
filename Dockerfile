FROM python:3.8-slim
WORKDIR /home/biolib
COPY requirements.txt .
RUN apt-get update && apt-get -y install dssp build-essential wget libjson-c-dev libjson-c3 libxml2 libxml2-dev libxml2-utils autoconf automake libtool
WORKDIR /opt
RUN wget http://freesasa.github.io/freesasa-2.0.3.tar.gz
RUN tar -xzf freesasa-2.0.3.tar.gz
WORKDIR /opt/freesasa-2.0.3
RUN ./configure CFLAGS="-fPIC -O2" --prefix=`pwd` --disable-xml
RUN make
RUN make install
WORKDIR /opt
RUN rm freesasa-2.0.3.tar.gz
RUN pip install -r requirements.txt
COPY predict.py .
COPY data/test.zip data/
ENTRYPOINT ["python3", "predict.py"]
