FROM rocker/tidyverse:4.3.2

# Install from CRAN
## Utility
RUN R -e "install.packages('parallel',dependencies=TRUE, repos='http://cran.rstudio.com/')"
## Dataframe
RUN R -e "install.packages('DT',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('data.table',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('dtplyr',dependencies=TRUE, repos='http://cran.rstudio.com/')"
## ggplot
RUN R -e "install.packages('viridis',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('igraph',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('ggupset',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('patchwork',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('ggstar',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('ggnetwork',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('gggenes',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('scales',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('hexbin',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('ggpointdensity',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('ggrepel',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('scatterpie',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('ggnewscale',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('png',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('svglite',dependencies=TRUE, repos='http://cran.rstudio.com/')"
## plotly
RUN R -e "install.packages('plotly',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('htmlwidgets',dependencies=TRUE, repos='http://cran.rstudio.com/')"
## tree
RUN R -e "install.packages('ape',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('tidytree',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('vegan',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('philentropy',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('pals',dependencies=TRUE, repos='http://cran.rstudio.com/')"
## rentrez
RUN R -e "install.packages('rentrez',dependencies=TRUE, repos='http://cran.rstudio.com/')"

# Install by BioManager
RUN R -e "install.packages('BiocManager',dependencies=TRUE, repos='http://cran.rstudio.com/')"
## Biostrings
RUN R -e "BiocManager::install('Biostrings')"
RUN R -e "BiocManager::install('coRdon')"
## ggtree
RUN R -e "BiocManager::install('treeio')"
RUN chmod -R a+w /usr/local/lib/R/site-library/ /usr/local/lib/R/library/ /usr/local/lib/R/doc/html/
# ggtree requires ggplot2 > 3.3.6 (rocker/tidyverse includes ggtree 3.3.6)
#RUN R -e "install.packages('tidyverse',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install('ggtree')"
RUN R -e "BiocManager::install('ggtreeExtra')"
## ggmsa
RUN apt-get update && apt-get install -y wget libxt6
RUN apt-get update && apt-get install -y libxml2-dev
RUN apt-get update && apt install -y libpng-dev
RUN apt-get update && apt install -y libfontconfig1-dev
RUN apt-get update && apt install -y libproj-dev
RUN R -e "BiocManager::install('ggmsa')"
RUN R -e "BiocManager::install('msa')"
## plyranges
RUN apt-get update && apt-get install -y libbz2-dev
RUN R -e "BiocManager::install('plyranges')"
## bio3d
RUN R -e "install.packages('bio3d',dependencies=TRUE, repos='http://cran.rstudio.com/')"

# Install by devtools
#RUN R -e "install.packages('usethis',dependencies=TRUE, repos='http://cran.rstudio.com/')"
#RUN R -e "install.packages('devtools',dependencies=TRUE, repos='http://cran.rstudio.com/')"
