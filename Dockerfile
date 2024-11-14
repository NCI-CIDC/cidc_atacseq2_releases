# Base image 
FROM amd64/ubuntu:20.04
MAINTAINER Frankie Parks <frankie.parks@bioappdev.org>

## Run all root level commands 1st
# Create a user and group used to launch processes
RUN groupadd -r pipeline -g 1000 && useradd -u 1000 -r -g pipeline -m -d /home/pipeline -s /bin/bash -c "pipeline user" pipeline && chmod 755 /home/pipeline && chown pipeline:pipeline /home/pipeline

#Create mount point directory
RUN mkdir -p /media/analysis && chmod 777 /media/analysis

## Update to latest OS packages
RUN apt-get autoclean && apt-get update && apt-get -y dist-upgrade
RUN apt-get -y install wget vim curl bc less python python3 python3-pip


USER pipeline
WORKDIR /home/pipeline
RUN pwd

RUN python -V

## Download and install gcloud utils
RUN curl -k -sSL https://sdk.cloud.google.com | bash
ENV PATH $PATH:/home/pipeline/google-cloud-sdk/bin
RUN gcloud components update --quiet

## Download and install miniconda
#RUN wget "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh" && bash Mambaforge-$(uname)-$(uname -m).sh -b -u
RUN wget "https://github.com/conda-forge/miniforge/releases/download/23.11.0-0/Mambaforge-$(uname)-$(uname -m).sh" && bash Mambaforge-$(uname)-$(uname -m).sh -b -u

# export miniconda path
ENV PATH="/home/pipeline/mambaforge/bin:$PATH"


# Set the working directory for pipeline code
RUN mkdir cidc_atac

WORKDIR /home/pipeline/cidc_atac

#Copy in all of the files/directories/etc.. from the project into the container
COPY . .
USER root
RUN chown -R pipeline:pipeline /home/pipeline/cidc_atac
USER pipeline
#ENV PATH="/home/pipeline/mambaforge/bin:$PATH"


## Install and initialize a mamba environment
RUN mamba init

## Add conda channels and set priorities
RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge
RUN conda config --set channel_priority flexible
RUN conda env create -f=environment.yaml
#RUN conda create -n nci-cli -c bioconda python=3.11 

#RUN conda activate cidc_atac

#ENTRYPOINT ["conda", "activate", "cidc_atac"]

CMD ["/bin/bash"]
