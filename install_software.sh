#!/bin/bash
#installs needed software
sudo sh -c 'echo "deb https://cran.cnr.berkeley.edu/bin/linux/ubuntu trusty/" >> /etc/apt/sources.list';
sudo apt-get update && sudo apt-get dist-upgrade -y;
sudo apt-get install build-essential libboost-all-dev git zlib1g-dev libcairo2-dev libxt-dev libx11-dev libncurses5-dev libcurl3 php5-curl r-base-dev bedtools samtools bwa -y python-numpy cython default-jre ant --force-yes;
mkdir ucsc ;
cd ucsc;
wget -m -e robots=off --no-parent http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/;
cd hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64;
rm -rf hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/index* hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat;
for f in hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/*; do 
 chmod +x $f; 
 sudo cp $f /usr/bin; 
done;
cd ..
git clone https://github.com/taoliu/MACS;
cd MACS;
sudo python setup_w_cython.py install;
cd bin;
chmod +x *;
sudo cp * /usr/bin;
cd ../..;
/usr/bin/Rscript install.R;
exit 0
