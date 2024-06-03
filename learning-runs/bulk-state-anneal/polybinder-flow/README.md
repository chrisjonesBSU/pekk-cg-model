# PolyBinder Flow 

This repository contains the scripts needed in order to use [Signac and Signac-flow](https://docs.signac.io/en/latest/) to create a state point space, and initialize simulations on a cluster.

## Instructions:

### Install and follow the set up instructions for PolyBinder:

[PolyBinder](https://github.com/cmelab/polybinder)

### Clone this repo:

`git clone git@github.com:cmelab/polybinder-flow.git`

### Set up the parameter space:  
- Edit the `init.py` file in the `src` directory with the desired statepoints.
- Run `python src/init.py` to generate the signac workspace

## Submit simulations:
- Run `python src/project.pyv submit` to run plolybinder (generate systems and initialize simulations)


