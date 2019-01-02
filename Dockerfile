# basic
FROM python:3.7.0

# maintainer information
MAINTAINER ygidtu ygidtu@gmail.com

# copy file to image
COPY sashimi /


RUN pip install pip -U
RUN pip config set global.index-url https://pypi.tuna.tsinghua.edu.cn/simple

# install requirements
RUN pip install -r requirements.txt

# Run the container as an executable
ENTRYPOINT ["python3", "main.py"]
