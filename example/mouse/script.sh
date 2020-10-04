#!/bin/bash
##
## cf https://software.intel.com/en-us/articles/installing-intel-free-libs-and-python-apt-repo

cd /tmp

#wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
# now install that key
sudo apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
# now remove the public key file exit the root shell
rm GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
exit
## all products:
#sudo wget https://apt.repos.intel.com/setup/intelproducts.list -O /etc/apt/sources.list.d/intelproducts.list
## just MKL
sudo sh -c 'echo deb https://apt.repos.intel.com/mkl all main > /etc/apt/sources.list.d/intel-mkl.list'

sudo apt-get update
sudo apt-get install intel-mkl-64bit-2019.5-075   ## wants 500+mb :-/  installs to 1.8 gb :-/

## update alternatives
update-alternatives --install /usr/lib/x86_64-linux-gnu/libblas.so     libblas.so-x86_64-linux-gnu      /opt/intel/mkl/lib/intel64/libmkl_rt.so 50
update-alternatives --install /usr/lib/x86_64-linux-gnu/libblas.so.3   libblas.so.3-x86_64-linux-gnu    /opt/intel/mkl/lib/intel64/libmkl_rt.so 50
update-alternatives --install /usr/lib/x86_64-linux-gnu/liblapack.so   liblapack.so-x86_64-linux-gnu    /opt/intel/mkl/lib/intel64/libmkl_rt.so 50
update-alternatives --install /usr/lib/x86_64-linux-gnu/liblapack.so.3 liblapack.so.3-x86_64-linux-gnu  /opt/intel/mkl/lib/intel64/libmkl_rt.so 50

echo "/opt/intel/lib/intel64"     >  /etc/ld.so.conf.d/mkl.conf
echo "/opt/intel/mkl/lib/intel64" >> /etc/ld.so.conf.d/mkl.conf
ldconfig

echo "MKL_THREADING_LAYER=GNU" >> /etc/environment
