### Base
FROM jupyter/datascience-notebook:aarch64-python-3.10.5

### apt-get installs
USER root
RUN sudo apt-get update && \
    apt-get install -y vim tree wget unzip bash-completion && \
    apt-get install -y libudunits2-dev && \
    apt-get install -y libnlopt-dev && \
    apt-get install -y default-jre less

### mamba install some R dependencies (needed for ggseg, ggpubr, images)
RUN mamba install -y libgdal nlopt r-units r-rcpp r-sf rpy2

### R packages
COPY r-reqs.R ./
RUN Rscript r-reqs.R


### mamba install 'vtk' dependency for python
RUN mamba install -y vtk

### Python reqs
COPY python-reqs.txt ./
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
ENV PATH $HOME/.cargo/bin:$PATH
RUN pip3 install gseapy
RUN pip3 install -r python-reqs.txt


### ENIGMA
#RUN git clone https://github.com/MICA-MNI/ENIGMA.git
#RUN python3 ENIGMA/setup.py install
#ENV PYTHONPATH "${PYTHONPATH}:/home/jovyan/ENIGMA"


### Connectome Workbench
RUN wget -O- http://neuro.debian.net/lists/focal.us-tn.libre | sudo tee /etc/apt/sources.list.d/neurodebian.sources.list
RUN sudo apt-get install -y gnupg
RUN sudo apt-key adv --recv-keys --keyserver keyserver.ubuntu.com 0xA5D32F012649A5A9
RUN sudo apt install ca-certificates
RUN sudo apt-get update
RUN sudo apt-get install -y connectome-workbench


ENV MPLCONFIGDIR=~/tmp/matplotlib

### [Legacy – switched to VS Code]
### Jupyter lab and extensions
#RUN pip3 install --upgrade jupyterlab
#RUN pip3 install -U "nbclassic>=0.2.8"
#RUN pip3 install lckr-jupyterlab-variableinspector
#RUN pip3 install aquirdturtle_collapsible_headings
#RUN jupyter labextension install @hokyjack/jupyterlab-monokai-plus 
##RUN pip3 install jupyterlab-lsp
#RUN pip3 install 'python-lsp-server[all]'
#RUN pip3 install jupyterlab_sublime
#RUN jupyter labextension install @techrah/text-shortcuts --no-build
#RUN jupyter labextension install @jupyterlab/debugger --no-build
#RUN jupyter lab build --dev-build=False --minimize=False
