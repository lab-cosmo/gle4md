gle4md: Generalized Langevin Equations for Molecular Dynamics
=============================================================


Utilities to generate input parameters for Generalized Langevin Equation thermostats

This repository collects tools to fit a set of drift and diffusion matrices for a generalized
Langevin equation, typically to be used in the context of molecular dynamics to enhance sampling
efficiency, or to enforce a non-equilibrium, frequency-dependent energy distribution. 

The code is fully functional and has been used in many papers, but is not exactly well 
documented and easy to use - and unlikely to get cleaned up and documented in the near future.

If you are interested in using it, and/or to contribute to the process of making it more
accessible, please get in touch with the authors!

Compilation
-----------

In order to compile gle4md you need to first download and compile the `toolbox` library
(you may have to adjust the settings in the `make.in` file in there

```
git clone git@github.com:lab-cosmo/toolbox.git
cd toolbox/src
make
```

then you should copy `make.in.example` to `make.in`, and tweak it to your system 
(e.g. setting the location where you have compiled the `toolbox` library. 
Then you should be able to just

```
cd src
make
```
