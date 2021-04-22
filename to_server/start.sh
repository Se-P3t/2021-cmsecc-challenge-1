#!/bin/bash

CURDIR=`pwd`

apt-get update && apt-get upgrade -y && apt autoremove -y
rm /usr/local/bin/python
ln -s `which python3` /usr/local/bin/python
apt-get install python3-pip -y
mkdir -p ~/.cache/pip/http
python -m pip install --upgrade pip
chmod a+w $CURDIR

apt-get install screen htop

apt-get install autoconf libtool m4 libgmp3-dev libmpfr-dev libqd-dev pkg-config -y
cd $CURDIR && git clone https://github.com/fplll/fplll.git
cd $CURDIR/fplll && ./autogen.sh && ./configure --with-max-parallel-enum-dim=160 && make -j3 && make install && ldconfig

cd $CURDIR && git clone https://github.com/fplll/fpylll.git
python -m pip install Cython
python -m pip install -r $CURDIR/fpylll/requirements.txt
python -m pip install -r $CURDIR/fpylll/suggestions.txt
cd $CURDIR/fpylll && python setup.py build_ext -j 6
cd $CURDIR/fpylll && python setup.py install && ldconfig

cd $CURDIR && git clone https://github.com/fplll/g6k.git
python -m pip install -r $CURDIR/g6k/requirements.txt
sed -i "s/maxsievingdim=128/maxsievingdim=1024/g" $CURDIR/g6k/rebuild.sh
cd $CURDIR/g6k && ./rebuild.sh -f && python setup.py install && ldconfig
