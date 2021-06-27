#/bin/bash
set -e

echo "use python3"
{
    rm -f /usr/local/bin/python
    ln -s `which python3` /usr/local/bin/python
} &> /dev/null


echo "update"
{
    apt-get update && apt-get upgrade -y
    apt-get install screen gdb -y

    python -m pip install --upgrade pip
} &> /dev/null


echo "install fplll"
{
    apt-get install autoconf libtool m4 libgmp3-dev libmpfr-dev libqd-dev pkg-config -y

    if [ -f /content/drive/MyDrive/fplll.tar.gz ]; then
        cd /root && tar -xpz -f /content/drive/MyDrive/fplll.tar.gz -C .
    else
        cd /root && git clone https://github.com/fplll/fplll.git
        cd /root/fplll && ./autogen.sh && ./configure && make -j2
        cd /root && tar -cpz -f /content/drive/MyDrive/fplll.tar.gz fplll/
    fi

    cd /root/fplll && make install && ldconfig
} &> /dev/null


echo "install fpylll"
{
    python -m pip install Cython

    if [ -f /content/drive/MyDrive/fpylll.tar.gz ]; then
        cd /root && tar -xpz -f /content/drive/MyDrive/fpylll.tar.gz -C .
        python -m pip install -r /root/fpylll/requirements.txt
        python -m pip install -r /root/fpylll/suggestions.txt
    else
        cd /root && git clone https://github.com/fplll/fpylll.git
        python -m pip install -r /root/fpylll/requirements.txt
        python -m pip install -r /root/fpylll/suggestions.txt
        cd /root/fpylll && python setup.py build_ext
        cd /root && tar -cpz -f /content/drive/MyDrive/fpylll.tar.gz fpylll/
    fi

    cd /root/fpylll && python setup.py install && ldconfig
} &> /dev/null


echo "install G6K-GPU-Tensor"
{
    if [ -f /content/drive/MyDrive/G6K-GPU.tar.gz ]; then
        #rm -rf /root/G6K-GPU/
        cd /root && tar -xpz -f /content/drive/MyDrive/G6K-GPU.tar.gz -C .
    else
        cd /root && git clone https://github.com/Se-P3t/G6K-GPU.git
        cd /root && tar -cpz -f /content/drive/MyDrive/G6K-GPU.tar.gz G6K-GPU/
    fi

    python -m pip install -r /root/G6K-GPU/requirements.txt
    cd /root/G6K-GPU && ./rebuild.sh -f -y && python setup.py install && ldconfig
} &> /dev/null


echo "logfile /tmp/screen.log" > /root/.screenrc
echo "DONE"
