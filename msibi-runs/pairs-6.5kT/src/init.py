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

def get_parameters():
    '''Use the listed parameters below to set up
    your MSIBI instructions.

    '''
    parameters = OrderedDict()
    parameters["head_correction"] = ["linear"]
    parameters["smooth"] = [True] # Whether or not to smooth the distributions
    parameters["smooth_potential"] = [True]
    parameters["smoothing_window"] = [5, 7, 9]
    parameters["integrator"] = ["hoomd.md.integrate.nvt"] # Hoomd integrator type
    parameters["integrator_kwargs"] = [{"tau": 0.01}] # dict of integrator kwargs
    parameters["nlist"] = ["hoomd.md.nlist.cell"]
    parameters["nlist_exclusions"] = [["1-2", "1-3"]]
    parameters["dt"] = [0.0001]
    parameters["gsd_period"] = [5e3] # Num of steps between gsd snapshots
    parameters["iterations"] = [20] # Num of MSIBI iterations to perform
    parameters["state_alphas"] = [
            [1.0, 1.0, 1.0, 1.0],
    ]
    parameters["n_steps"] = [
            #5e5,
            1e6,
            #3e6,
            5e6,
    ]
    parameters["optimize"] = ["pairs"] # Choose which potential to optimize

    # These parameters below are only needed when optimizing pair potentials
    parameters["rdf_exclude_bonded"] = [True] # Exclude pairs on same molecule
    parameters["r_switch"] = [None] # Distance value to apply tail correction

    # Add state points to use during MSIBI
    parameters["states"] = [
            # Evenly Weighted, with 1.0
        [
            {"name":"A",
            "kT":6.5,
            "target_trajectory":"two-chains.gsd",
            "max_frames": 150,
            "target_frames": 200,
            "exclude_bonded": False
            },

            {"name":"B",
            "kT":5.0,
            "target_trajectory":"5.0kT-1.27den.gsd",
            "max_frames": 100,
            "target_frames": 200,
            "exclude_bonded": True
            },

            {"name":"C",
            "kT":5.0,
            "target_trajectory":"5.0kT-1.38den.gsd",
            "max_frames": 100,
            "target_frames": 200,
            "alpha":1.0,
            "exclude_bonded": True
            },

            {"name":"D",
            "kT":6.5,
            "target_trajectory":"6.5kT-1.27den.gsd",
            "max_frames": 100,
            "target_frames": 200,
            "alpha":1.0,
            "exclude_bonded": True
            },
        ],

    ]

    # Add parameters needed to create Pair objects
    parameters["pairs"] = [
        [
            {"type1":"E",
             "type2":"E",
             "form": "table",
             "kwargs": {
                 "n_points": 101,
                 "epsilon": 1.5,
                 "sigma": 1,
                 "r_max": 4.0,
                 "r_min": 1e-3
              }
             },

            {"type1":"E",
             "type2":"K",
             "form": "table",
             "kwargs": {
                 "n_points": 101,
                 "epsilon": 1.5,
                 "sigma": 1,
                 "r_max": 4.0,
                 "r_min": 1e-3
              }
             },

            {"type1":"K",
             "type2":"K",
             "form": "table",
             "kwargs": {
                 "n_points": 101,
                 "epsilon": 1.5,
                 "sigma": 1,
                 "r_max": 4.0,
                 "r_min": 1e-3
              }
             },
        ]

	]

    # Add parameters needed to create Bond and Angle objects
    parameters["bonds"] = [
            [
                {"type1":"E",
                 "type2":"K",
                 "form": "harmonic",
                 "kwargs": {"k": 850, "l0": 1.47}
                 },

                {"type1":"K",
                 "type2":"K",
                 "form": "harmonic",
                 "kwargs": {"k": 1450, "l0": 1.53}
                 }
            ],
    ]

    parameters["angles"] = [
            [
                {"type1":"E",
                "type2":"K",
                "type3": "K",
                "form": "file",
                "kwargs": {"file_path": "E-K-K_angle_1.0.txt"}
                },

                {"type1":"K",
                "type2":"E",
                "type3": "K",
                "form": "file",
                "kwargs": {"file_path": "K-E-K_angle_1.0.txt"}
                },
            ]
    ]

    parameters["dihedrals"] = [
            [
                {"type1":"E",
                "type2":"K",
                "type3": "K",
                "type4": "E",
                "form": "harmonic",
                "kwargs": {"phi0": 0, "k": "16", "d": -1, "n":1}
                },

                {"type1":"K",
                "type2":"E",
                "type3": "K",
                "type4": "K",
                "form": "harmonic",
                "kwargs": {"phi0": 0, "k": "12", "d": -1, "n":1}
                },

            ]
    ]
    return list(parameters.keys()), list(product(*parameters.values()))

def main():
    project = signac.init_project("msibi")
    param_names, param_combinations = get_parameters()
    # Create the generate jobs
    for params in param_combinations:
        parent_statepoint = dict(zip(param_names, params))
        parent_job = project.open_job(parent_statepoint)
        parent_job.init()
        parent_job.doc.setdefault("integrator_kwargs", parent_job.sp.integrator_kwargs)
    project.write_statepoints()

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
