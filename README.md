# msmrd2:
MSM/RD: A general framework to couple Markov models of molecular kinetics with reaction-diffusion simulations. This is the second version of the code rewritten in c++ with python bindings and additional functionality. (Work in progress)

## Installation
```
git clone https://github.com/maojrs/msmrd2.git
cd msmrd2
git submodule update --init --recursive
python setup.py install
```

## To do:
- [ ] Extend overdampedLangevinMS and msmrdIntegrator to handle discrete MSMs (not required but could be a nice addition).
