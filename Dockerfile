FROM kubor/alpine-rdkit
RUN apk add --no-cache git
ENV PYTHONBUFFERED 1
RUN mkdir /code
WORKDIR /code
ADD requirements.txt /code/
RUN pip install -r requirements.txt
RUN pip install git+https://github.com/rvianello/razi.git
ADD . /code/
