FROM r-base:4.6.0

ENV CRANREPO="'https://pbil.univ-lyon1.fr/CRAN'"
ENV RLIBPATH="'/usr/local/lib/R/site-library/'"

# File Author / Maintainer
LABEL maintainer="alemeur@biophylo.com"

RUN apt-get update --fix-missing && apt-get install -y wget gcc bzip2 make zlib1g-dev libcurl4-gnutls-dev python3 curl pip git python-dev-is-python3
RUN pip3 install cyvcf2 numpy=="1.26.4" --break-system-packages

RUN R --silent -e "install.packages('BiocManager',repo=$CRANREPO, lib=$RLIBPATH)" \
	&& R --silent -e "install.packages('Rcpp',repo=$CRANREPO, lib=$RLIBPATH)" \
	&& R -e "install.packages('seqinr',repo=$CRANREPO, lib=$RLIBPATH)" \
	&& R --silent -e "install.packages('R6',repo=$CRANREPO, lib=$RLIBPATH)" \
	&& R --silent -e "install.packages('ape',repo=$CRANREPO, lib=$RLIBPATH)" \
	&& R --silent -e "BiocManager::install('Rhtslib',repo=$CRANREPO, lib=$RLIBPATH)" \
	&& R --silent -e "install.packages('optparse',repo=$CRANREPO, lib=$RLIBPATH)" \
	&& R --silent -e "BiocManager::install('Rcurl',repo=$CRANREPO, lib=$RLIBPATH)" \
	&& R --silent -e "BiocManager::install('Biostrings',repo=$CRANREPO, lib=$RLIBPATH)" \
	&& R --silent -e "install.packages('dplyr',repo=$CRANREPO, lib=$RLIBPATH)" \
	&& R --silent -e "install.packages('maps',repo=$CRANREPO, lib=$RLIBPATH)" \
	&& R --silent -e "install.packages('jackalope',repo=$CRANREPO, lib=$RLIBPATH)"

RUN git clone "https://github.com/adrienlemeur/maketube.git" "/usr/local/bin/maketube" 	&& chmod +x -R "/usr/local/bin/maketube"
RUN export PATH="/usr/local/bin/maketube:$PATH"

ENTRYPOINT ["Rscript", "/usr/local/bin/maketube/maketube.R"]
