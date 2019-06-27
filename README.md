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

- [x] Finish working version of generalized MSM/RD integrator with arbitrary relative position/orientation.

- [ ] Implement first test case for simplest two dimer/two bound states model.

- [ ] Implement test case for 'patchy proteins'.

- [ ] Implement case for trimer and pentamer (requires additional modification of msmrdIntegrator, perhaps make a child class). 

- [ ] Review untested parts of code and write tests (not urgent since most complex and relevant parts already have extensive tests).

- [ ] Extend overdampedLangevinMS and msmrdIntegrator to handle discrete MSMs (not required but could be a nice addition).
