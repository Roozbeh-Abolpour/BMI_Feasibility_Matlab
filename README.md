# MATLAB Codes for BMI Feasibility Algorithms

This repository contains MATLAB implementations of several algorithms for bounded bilinear matrix inequality (BMI) feasibility problems. The code is intended for numerical comparison of BMI-oriented feasibility methods.

## Repository structure

```text
Bab/        Branch-and-bound implementation
Bender/     Benders-decomposition-based implementation
general/    Shared utility functions
ILMI/       Iterative LMI implementation
Mine/       Proposed separating-hyperplane/inner-loop implementation
Path/       Path-following implementation
```

## Main routines

Each main solver uses the interface

```matlab
[x,y,T] = solver_name(Ls,Bs,bounds,params)
```

where `solver_name` can be one of:

- `bab` for branch-and-bound;
- `bender` for Benders decomposition;
- `ilmi` for iterative LMI;
- `mine` for the proposed separating-hyperplane refinement method;
- `bpath` for the path-following method.

## Problem format

The input cell arrays `Ls` and `Bs` define matrix-valued constraints as MATLAB function handles. Typically,

```matlab
Ls{i} = @(x,y) ... ;
Bs{i} = @(x,y) ... ;
```

The vector variables are bounded by

```matlab
bounds.lbx = ...;
bounds.ubx = ...;
bounds.lby = ...;
bounds.uby = ...;
```

The structure `params` contains method-dependent stopping tolerances and time limits, for example:

```matlab
params.conv_ther      = 1e-4;
params.accept_ther    = 1e-6;
params.tmax           = 60;
params.init_violation = 1e3;
params.sep_cap        = -1e-4;
params.path_radius    = 1e-1;
params.conv_length    = 1e-3;
```

Not every field is used by every method.

## Dependencies

The implementation uses:

- MATLAB;
- YALMIP;
- MOSEK as the SDP solver.

Make sure the required folders are added to the MATLAB path before running the methods:

```matlab
addpath(genpath(pwd));
```

## Notes

These codes are research implementations and are provided for reproducibility of numerical experiments. Before using them in a new benchmark, check that the tolerances, solver settings, and time limits match the intended experimental protocol.
