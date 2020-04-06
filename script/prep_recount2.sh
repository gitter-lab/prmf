#!/bin/sh
set -e
cd ~/data
wget -O recount2_PLIER_data.zip https://ndownloader.figshare.com/files/10881866
unzip recount2_PLIER_data.zip
#docker run -it -v `pwd`:`pwd` -w `pwd` aabaker99/cogaps:latest Rscript prep_recount2.R
