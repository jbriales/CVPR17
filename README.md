# Code for Convex Global 3D Registration with Lagrangian Duality (CVPR 17)

This repository contains the implementation of the methods compared in

@inproceedings{briales17CVPR,
title = {{Convex Global 3D Registration with Lagrangian Duality}},
author = {Briales, Jesus and Gonzalez-Jimenez, Javier},
booktitle = {International Conference on Computer Vision and Pattern Recognition},
month = {jul},
year = {2017}
}

In this work we proposed a novel convex relaxation for registration of points to points, lines and planes.
Empirically, the relaxation results always tight and certifies global optimality.
Indeed, the framework is able to deal with any optimization problem that has a quadratic objective
on rotation matrix elements.


## Getting started

### Clone the repository
This repository include some dependencies as submodules,
so clone it with the `--recursive` option:
```
  git clone --recursive https://github.com/jbriales/CVPR17.git
```
If you already cloned it, you can still set the submodules with
```
  git submodule update --init --recursive
```

### Install the library
To use the provided code and methods, just run the `setup.m` script.
Note you should have installed *CVX* (available [here](http://cvxr.com/cvx/)) in the path:

- [Download CVX](http://cvxr.com/cvx/download/) for your platform
- Install CVX: Run `cvx_setup.m` from Matlab

### Run the examples
For a working example, see *example.m*.

