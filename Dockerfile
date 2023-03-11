FROM ubuntu

RUN apt update
RUN apt install python3-pip -y
RUN apt-get install bedtools
# RUN apt-get install libbz2-dev liblzma-dev
RUN pip3 install pandas
RUN pip3 install pybedtools
RUN pip3 install seaborn
RUN pip3 install matplotlib
RUN pip3 install swifter
RUN pip install matplotlib-venn
RUN pip3 install Flask
RUN pip3 install waitress

WORKDIR /app

COPY . .

CMD ["python3","waitress_server.py"]