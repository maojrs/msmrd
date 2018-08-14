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

- [x] Test sample overdamped Lanegvin integrator with pybindings.

- [x] Write header and source file for potentials.

- [x] Rewrite module for MSM and CTMSM in c++.

- [x] Rewrite previous code for diffusion with rotation and Markovian switch. 

- [ ] Write tests for randgen, msm, ctmsm, odLangevin and odLangevinMarkovSwitch classes

- [ ] Choose toy system with orientation dependent potential
