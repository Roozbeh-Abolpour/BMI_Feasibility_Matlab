# MATLAB Codes for BMI Feasibility Algorithms

This repository contains MATLAB implementations of several algorithms for bounded bilinear matrix inequality (BMI) feasibility problems. The codes are intended for reproducing and comparing the numerical performance of BMI-oriented feasibility methods.

## Repository structure

```text
Bab/        Branch-and-bound implementation
Bender/     Benders-decomposition-based implementation
general/    Shared utility functions
ILMI/       Iterative LMI implementation
Mine/       Proposed separating-hyperplane refinement implementation
Path/       Path-following implementation
PSDP/       Penalized semidefinite programming implementation
SDP/        Semidefinite programming relaxation
SOCP/       Second-order cone programming relaxation
```

## Main routines

Each main solver follows the interface

```matlab
[x,y,T] = solver_name(Ls,Bs,bounds,params)
```

where:

* `x` and `y` are the computed decision vectors;
* `T` is the elapsed solution time measured by the implementation;
* `solver_name` can be one of the following:

  * `bab` for branch-and-bound;
  * `bender` for Benders decomposition;
  * `ilmi` for iterative LMI;
  * `mine` for the proposed separating-hyperplane refinement method;
  * `bpath` for the path-following method;
  * `bsdp` for the semidefinite programming relaxation;
  * `bsocp` for the second-order cone programming relaxation;
  * `psdp` for the penalized semidefinite programming method.

## Problem format

The methods seek bounded vectors `x` and `y` satisfying the supplied linear matrix inequality and bilinear matrix inequality constraints.

The input cell arrays `Ls` and `Bs` define matrix-valued constraints as MATLAB function handles. Typically,

```matlab
Ls{i} = @(x,y) ... ;
Bs{i} = @(x,y) ... ;
```

The decision-variable bounds are specified by

```matlab
bounds.lbx = ...;
bounds.ubx = ...;
bounds.lby = ...;
bounds.uby = ...;
```

The structure `params` contains method-dependent stopping tolerances, algorithmic parameters, and time limits. For example,

```matlab
params.conv_ther      = 1e-4;
params.accept_ther    = 1e-6;
params.tmax           = 60;
params.init_violation = 1e3;
params.sep_cap        = -1e-4;
params.path_radius    = 1e-1;
params.conv_length    = 1e-3;
params.penalty_gain   = 1e2;
```

Not every field is used by every method. In particular, `penalty_gain` is used by the penalized SDP implementation.

## Dependencies

The implementations require:

* MATLAB;
* YALMIP;
* MOSEK as the semidefinite and conic optimization solver.

Make sure that YALMIP and MOSEK are installed and configured. Add the repository and all its subfolders to the MATLAB path before running the methods:

```matlab
addpath(genpath(pwd));
```

## Basic usage

After defining `Ls`, `Bs`, `bounds`, and the required fields of `params`, call a solver directly. For example,

```matlab
[x,y,T] = mine(Ls,Bs,bounds,params);
```

For the SDP relaxation,

```matlab
[x,y,T] = bsdp(Ls,Bs,bounds,params);
```

For the SOCP relaxation,

```matlab
[x,y,T] = bsocp(Ls,Bs,bounds,params);
```

For the penalized SDP method,

```matlab
[x,y,T] = psdp(Ls,Bs,bounds,params);
```

## Notes

These codes are research implementations provided for reproducibility of numerical experiments. Before using them in a new benchmark, verify that the stopping tolerances, solver settings, penalty parameters, and time limits match the intended experimental protocol.
