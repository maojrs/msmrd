# msmrd:
MSM/RD: A general framework to couple Markov models of molecular kinetics with reaction-diffusion simulations. This is the second version of the code (v2.0.0) rewritten in c++ with python bindings and additional functionality. The first version (v1.0.0) was a much simpler implementation written in python, and it can be found [here](https://github.com/markovmodel/msmrd). (Work in progress)

## Installation
```
git clone https://github.com/maojrs/msmrd2.git
cd msmrd2
git submodule update --init --recursive
python setup.py install
```

We recommend creating a conda environment for msmrd

## Software dependencies 
Check the environment.yml for the specific dependencies in the conda environment. In general, the following dependencies are required:
- GCC 10.2.0
- cmake 3.19.4
- HDF5 1.10.4
- python 3.8.3
- numpy 1.18.1
- scipy 1.4.1
- matplotlib 3.2.1
- h5py 2.10.0

Versions below may work; versions above will most likely work. for the notebooks, we recommend also installing [Jupyter notebooks](https://jupyter.org/)
