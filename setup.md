## env setup ##

The scripts require the SageMath environment to run, it's helpful to start with a Conda env. We prefer to choose [mamba](https://github.com/mamba-org/mamba) as the package manager.

```sh
mamba create -y -n g6k-env sage python=3.9
```

Sage alreadly have FPyLLL inside. To get the sieve work, we only need to install G6K.

```sh
git clone https://github.com/fplll/g6k
cd g6k
mamba install -y -n g6k-env --file requirements.txt

conda activate g6k-env
python setup.py clean
./rebuild.sh
python setup.py install
#python -m pytest
conda deactivate
```

### [HPLLL](https://perso.ens-lyon.fr/gilles.villard/hplll) ###

```sh
./configure --enable-omp --with-gmp=$CONDA_PREFIX --with-fplll=$CONDA_PREFIX --prefix=$CONDA_PREFIX

g++ -O3 -Wall -pthread -fopenmp hlll.cpp -o hlll -lmpfr -lgmp -lquadmath -I$CONDA_PREFIX/include -L$CONDA_PREFIX/lib
```
