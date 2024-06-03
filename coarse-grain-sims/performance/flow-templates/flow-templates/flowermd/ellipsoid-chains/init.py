#!/usr/bin/env python
"""Initialize the project's data space.

Iterates over all defined state points and initializes
the associated job workspace directories.
The result of running this file is the creation of a signac workspace:
    - signac.rc file containing the project name
    - signac_statepoints.json summary for the entire workspace
    - workspace/ directory that contains a sub-directory of every individual statepoint
    - signac_statepoints.json within each individual statepoint sub-directory.

"""

import signac
import flow
import logging
from collections import OrderedDict
from itertools import product


def get_parameters():
    ''''''
    parameters = OrderedDict()
    # Beads and system:
    parameters["num_mols"] = [[1]]
    parameters["lengths"] = [[25]]
    parameters["density"] = [0.01] # number density
    parameters["bead_mass"] = [1]

    # GB Params
    parameters["lpar"] = [1] # nm
    parameters["lperp"] = [0.5] # nm
    parameters["epsilon"] = [1.0] # kJ/mol

    # Bonds
    parameters["r0"] = [0.01]
    parameters["bond_k"] = [500]

    # Angles
    parameters["theta0"] = [2.8]
    parameters["theta_k"] = [200]

    # Define some simulation related parameters:
    parameters["kT"] = [0.5]
    parameters["n_steps"] = [1e7]
    parameters["shrink_kT"] = [4.0]
    parameters["shrink_n_steps"] = [1e6]
    parameters["shrink_period"] = [10000]
    parameters["r_cut"] = [2.5]
    parameters["dt"] = [0.0005]
    parameters["tau_kT"] = [100] # Used as a multiple of dt
    parameters["gsd_write_freq"] = [2e4]
    parameters["log_write_freq"] = [2e3]
    parameters["seed"] = [24]
    return list(parameters.keys()), list(product(*parameters.values()))


def main():
    project = signac.init_project() # Set the signac project name
    param_names, param_combinations = get_parameters()
    # Create the generate jobs
    for params in param_combinations:
        statepoint = dict(zip(param_names, params))
        job = project.open_job(statepoint)
        job.init()
        job.doc.setdefault("nvt_done", False)
        job.doc.setdefault("sample_done", False)


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
