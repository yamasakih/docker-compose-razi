FROM jupyter/scipy-notebook 

USER root
RUN conda install -c rdkit rdkit
ENV PYTHONBUFFERED 1

RUN mkdir /code
WORKDIR /code
ADD requirements.txt /code/
RUN pip install -r requirements.txt
RUN pip install git+https://github.com/rvianello/razi.git

WORKDIR /home/jovyan
USER jovyan

