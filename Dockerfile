FROM jupyter/scipy-notebook 

ENV PYTHONBUFFERED 1

USER root
RUN conda install -c rdkit rdkit
RUN conda install -c conda-forge molvs

RUN mkdir /code
WORKDIR /code
ADD requirements.txt /code/
RUN pip install -r requirements.txt
RUN pip install git+https://github.com/rvianello/razi.git

WORKDIR /home/jovyan
USER jovyan

