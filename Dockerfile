FROM python:3.8-slim
WORKDIR /home/biolib
COPY requirements.txt .
RUN apt-get update && apt-get -y install dssp build-essential wget libjson-c-dev libjson-c3 libxml2 libxml2-dev libxml2-utils autoconf automake libtool
WORKDIR /home/biolib
RUN pip install -r requirements.txt
COPY predict.py .
COPY featcomputers.py .
COPY first_model.bin data/
COPY lasso_model_v1.bin data/
COPY data/test.zip data/
ENTRYPOINT ["python3", "predict.py"]
