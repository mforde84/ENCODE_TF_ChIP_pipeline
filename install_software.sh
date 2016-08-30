#!/bin/bash
#installs needed software
sudo sh -c 'echo "deb https://cran.cnr.berkeley.edu/bin/linux/ubuntu trusty/" >> /etc/apt/sources.list'
sudo apt-get update && sudo apt-get dist-upgrade -y
sudo apt-get install build-essential libboost-all-dev git zlib1g-dev libcairo2-dev libxt-dev libx11-dev libncurses5-dev libcurl3 php5-curl r-base-dev bedtools samtools bwa -y python-numpy cython default-jre ant --force-yes
git clone https://github.com/taoliu/MACS
cd MACS
sudo python setup_w_cython.py install
cd bin
chmod +x *
sudo cp * /usr/bin
/usr/bin/Rscript install.R 
exit 0
