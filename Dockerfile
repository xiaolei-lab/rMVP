FROM r-base:3.6.3

RUN Rscript -e "update.packages(ask=F); install.packages('rMVP', repos='https://mirrors.tuna.tsinghua.edu.cn/CRAN/')"