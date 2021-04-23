FROM python:3.8-slim
WORKDIR /home/biolib
COPY requirements.txt .
RUN pip install -r requirements.txt
RUN apt-get install dssp
COPY predict.py .
ENTRYPOINT ["python3", "predict.py"]
