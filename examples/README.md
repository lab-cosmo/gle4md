GLE fitting examples
====================

This folder contains a couple of examples of GLE fitting inputs. 
The documentation is sadly non-existent, so here you can get at least
a starting point to experiment. 

Fitting the friction kernel spectrum K(w)
-----------------------------------------

`kw-fit.in` describes the target for a GLE fitting exercise in which we try
to match the Fourier transform of the memory kernel to a target spectral shape,
rather arbitrarily taken to be the velocity autocorrelation spectrum of liquid 
water. The target function is stored in `kw-tgt.dat`, and you can use `kw-gen.sh`
to subsample it and format it as a target definition that can be pasted in the 
input -- looking like this

```
{ frequency 7.34964e-08 values { kw = 2.49563e+00 1.0e+0 log 1  } }
{ frequency 5.87971e-07 values { kw = 2.71662e+00 1.0e+0 log 1  } }
{ frequency 1.98440e-06 values { kw = 3.02201e+00 1.0e+0 log 1  } }
```

You can use `gle-deltas.py` to prepare an initial guess for the A matrix that 
matches the target. Type `kw-tgt.dat` in the `Plot file` box, and adjust the sliders
until the K(w) curve resembles the target - you will need to have `gle-analyze`
in your path, and tick the `Auto-compute` box. 
Then you can copy and paste the content of `deltafit.a` into the `restart` field
in `kw-fit.in`. 

To run the fitting, you should do

```
$ gle-fit < kw-fit.in &> kw-fit.log
```

stare at the hypnotic flood of debug mess printed to the log, and then find the
optimized drift matrix in `kw-delta.A`. You can then compute the predicted GLE
properties by

```
$ gle-analyze -wi 1e-6 -wf 1e2 -a kw-delta.A > kw-fit.dat
```

and compare K(w) (which is column 5) with the target.


Fitting a PIGLET matrix
-----------------------

The input `piglet-fit.in` is an example of the fitting of a set of matrices for
path integrals + GLE simulations. You will see it's rather complicated, as we are 
trying to find a balance between many desiderata. A brief explanation of the 
target lines, that look like this

```
{ frequency 0.000200 values { cqq = 11.999993 1.0e+0 log 1   cpp = 11.999993 1.0e+0 log 1   kq2 = 1.000000 0.006072 log 1    cqqdt < 11.999993 1.0e+0 log 1  kw = 0.000200 0.006072 log 1 } }
```

goes as follows. `cpp` and `cqq` are the stationary fluctuations in <p2> and <q2>. 
These are set to match the PIGLET targets computed in the 
[PIGLET paper](https://doi.org/10.1103/PhysRevLett.109.100604).
Then we also ask for the sampling of q2 to be as efficient as possible, as indicated
by the target `kq2 = 1`. Then we ask for the dynamics not to be crazily sensitive on
the time step (as it can happen if it contains very high noise components) with
the `cqqdt` target (the time step corresponding to it is given at the beginning of 
the file, `finite-dt      1`). Finally, for good measure, we ask that `kw = omega`,
which tends to be a good balanced GLE requirement. Note that the weight of the 
different requirements is modulated by the global weights in `pweight`. 

So try to 
```
$ gle-fit < piglet-fit.in &> piglet-fit.log
```
and make some sense out of it. It's not a fantastic introduction to the art of GLE
fitting, but better than nothing.

