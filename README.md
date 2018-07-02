# msmrd2:
MSM/RD: A general framework to couple Markov models of molecular kinetics with reaction-diffusion simulations. This is the second version of the code rewritten in c++ with python bindings and additional functionality. (Work in progress)

## Installation
```
git clone https://github.com/markovmodel/msmrd2.git
cd msmrd2
git submodule update --init --recursive
python setup.py install
```

## To do:

- [ ] Test sample overdamped Lanegvin integrator with pybindings.

- [ ] Create header file for potentials.

- [ ] Rewrite module for MSM and CTMSM in c++.

- [ ] Rewrite previous code for diffusion with rotation and Markovian switch. 
