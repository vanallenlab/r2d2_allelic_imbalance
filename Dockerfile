FROM ubuntu:17.10

WORKDIR /

RUN apt-get update && \
	apt-get upgrade -y && \
	apt-get install python-tk -y && \
	apt-get install -y wget vim bzip2 less


RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda3-4.5.4-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    /opt/conda/bin/conda clean -tipsy && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

ENV PATH /opt/conda/bin:$PATH

RUN conda update conda -y

RUN pip install pandas
RUN pip install scipy

COPY r2d2v2.py /
COPY clinvar_variant_summary.lite.txt /
