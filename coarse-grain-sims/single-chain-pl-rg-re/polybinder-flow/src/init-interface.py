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
import logging
from collections import OrderedDict
from itertools import product
import numpy as np

def get_parameters():
    '''
    Parameters:
    -----------

    System generation parameters:
    -----------------------------


    Simulation parameters:
    ----------------------


    ------------
    Other Notes:
    ------------
    All temperatures are entered as reduced temperature units

    If you only want to run a quench simulation
        Comment out kT_anneal, anneal_sequence lines

    If you only want to run an anneal simulation
        Comment out kT_quench and n_steps lines

    Don't forget to change the name of the project
    project = signac.init_project("project-name")
    '''

    parameters = OrderedDict()
    # System generation parameters:
    parameters["signac_project"] = [
            "/home/erjank_project/chrisjones/tensile/make_slabs"
        ] # Path to projec that contains slabs.
                                          # Leave as None if they are in this current project
    parameters["signac_args"] = [["457240ab30858158b22ffceeda4bc813"] # A way for signac to find the slab .gsd file(s)
                                ]   # Can be a job ID or a dictionary of state points

    parameters["slab_file"] = [None]  # Full path to .gsd file(s)

    parameters["interface_gap"] = [0.1]
    parameters["weld_axis"] = ["z"]
    parameters["reference_distance"] = [0.339]
    parameters["forcefield"] = ['gaff']
    parameters["remove_hydrogens"] = [True]
    parameters["system_seed"] = [24]
    # Simulation parameters
    parameters["tau_kt"] = [0.1]
    parameters["tau_p"] = [None]
    parameters["pressure"] = [None]
    parameters["dt"] = [0.001]
    parameters["r_cut"] = [2.5]
    parameters["neighbor_list"] = ["cell"]
    parameters["sim_seed"] = [42]
    parameters["walls"] = [[0,0,1]]
    parameters["system_type"] = ["interface"] # Don't change this
    parameters["procedure"] = [#"quench",
                              "anneal"
                              ]

        # Quench related params:
    #parameters["kT_quench"] = [1.5] # Reduced Temp
    #parameters["n_steps"] = [1e7]

        # Anneal related params
    parameters["kT_anneal"] = [
                               [5.5, 2.0]
                              ] # List of [initial kT, final kT] Reduced Temps
    parameters["anneal_sequence"] = [
                                     [2e5, 1e5, 3e5, 5e5, 5e5, 1e5] # List of lists (n_steps)
                                    ]
    parameters["schedule"] = [None]
    return list(parameters.keys()), list(product(*parameters.values()))

custom_job_doc = {} # keys and values to be added to each job document created
                    # leave blank to create for job doc entries
def main():
    project = signac.init_project("project")
    param_names, param_combinations = get_parameters()
    # Create the generate jobs
    for params in param_combinations:
        parent_statepoint = dict(zip(param_names, params))
        parent_job = project.open_job(parent_statepoint)
        parent_job.init()
        try:
            parent_job.doc.setdefault("steps", parent_statepoint["n_steps"])
        except:
            parent_job.doc.setdefault("steps", np.sum(parent_statepoint["anneal_sequence"]))
            parent_job.doc.setdefault("step_sequence", parent_statepoint["anneal_sequence"])
        if parent_job.sp['signac_args'] is not [None]:
            parent_job.doc.setdefault("use_signac", True)
            parent_job.doc.setdefault("slab_files", parent_job.sp['signac_args'])
        elif parent_job.sp['slab_file'] is not [None]:
            parent_job.doc.setdefault("use_signac", False)
            parent_job.doc.setdefault("slab_files", parent_job.sp['slab_file'])

    if custom_job_doc:
        for key in custom_job_doc:
            parent_job.doc.setdefault(key, custom_job_doc[key])

    project.write_statepoints()

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
