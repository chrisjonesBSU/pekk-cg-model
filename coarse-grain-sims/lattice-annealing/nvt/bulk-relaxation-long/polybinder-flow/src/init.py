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
    '''
    Parameters:
    -----------

    System generation parameters:
    -----------------------------
    molecule : str
        Name of the molecule used to build the system.
        Must match one of the json files in uli-init/compounds
    para_weight : float; between 0 and 1
        The relative amount of para conformations in the system
        1 = All para, 0 = All meta
    density : float
        The density of the system in g/cm^3.
        PEEK and PEKK are both around 1.3 - 1.4 g/cm^3
    n_compounds : list
        A list of the number of molecules of a given length
        Must be the same legnth as polymer_lengths list(s)
        Corresponds to the number of specific molecules
        at the same index position in polymer_lengths
        See pdi parameter
    polymer_lengths : list
        A list of the number of monomer units in a single molecule
        Must be the same legnth as n_compounds list(s)
        See pdi parameter
    sample_pdi : bool
        Instruct uli-init to generate a distribution using a combination
        of pdi, Mn, Mw. This will override n_compound and polymer_length
        parameters
    pdi : float
        A PDI (poly-dispersity index) value of the generated system.
        pdi = Mn/Mw
    Mn : int
        The most frequent polymer length of a polydisperse system
        Used in conjunction with pdi to determine distribution
        of polymer lengths in the system
    Mw : int
        The weight average of the polymer distribution.
    forcefield : str options are 'gaff' or 'opls'
        The forcefield type to use when calling Foyer
    mass_dist : str
        Specify the distribution to be used when sampling from a pdi
        Options are: 'weibull' or 'gaussian'
    walls : np.array like
        Use walls to create flat surfaces along one of the volume axes
        Walls on x-axis: [1,0,0]
        Walls on y-axis: [0,1,0]
        Walls on z-axis: [0,0,1]
        Leaving as None will result in PBC across all volume faces

    Simulation parameters:
    ----------------------


    ------------
    Other Notes:
    ------------
    All temperatures are entered as reduced temperature units

    If you want to sample from a PDI:
        Change the polymer length lines to [None]
        Change the n_compounds lines to [None]

    If you only want to run a quench simulation
        Comment out kT_anneal, anneal_sequence lines

    If you only want to run an anneal simulation
        Comment out kT_quench and n_steps lines

    Don't forget to change the name of the project

    '''
    parameters = OrderedDict()

    ### SYSTEM GENERATION PARAMETERS ###
    parameters["system_type"] = ["crystal"]
    parameters["molecule"] = [
            #'PEEK',
            'PEKK',
            #"PPS"
    ]
    parameters["para_weight"] = [1.0]
    parameters["ekk_weight"] = [
            1.0,
            0.8,
            0.7,
            0.6
    ]
    parameters["kek_weight"] = [1.0]
    parameters["bonds"] = [
            {
                "E-K": {"k": 850, "r0": 1.47},
                "K-K": {"k": 1450, "r0": 1.53},
            }
    ]
    parameters["dihedrals"] = [
            {
                "E-K-K-E": {"k": 16, "phi0": 0, "d": -1, "n": 1 },
                "K-E-K-K": {"k": 12, "phi0": 0, "d": -1, "n": 1 },
            }
    ]
    parameters["monomer_sequence"] = [None]
    parameters["density"] = [
            1.35,
    ]
    parameters["n_compounds"] = [[200]]
    parameters["polymer_lengths"] = [[20]]
    parameters["pdi"] = [None]
    parameters["Mn"] = [None]
    parameters["Mw"] = [None]
    parameters['mass_dist'] = ['weibull']
    parameters["charges"] = [None]
    parameters["forcefield"] = [None]
    parameters["remove_hydrogens"] = [True]
    parameters["system_seed"] = [42]
    parameters["box_constraints"] = [
            {"x": 0.606 * 11, "y": 0.767 * 11, "z": None},
    ]
    parameters["kwargs"] = [
            {"n": 10, "a": 0.606, "b": 0.767}
    ]

    ### SIM FROM RESTART PARAMETERS ###

	# Path to the signac project to use
    parameters["signac_project"] = [None]
	# A way for signac to find the specific state point to use
	# Can be a job ID or a dictionary of a state point
    parameters["signac_args"] = [None]
	# Give the full path to the restart.gsd file instead of using signac
    parameters["restart_file"] = [None]

    ### COARSE-GRAINING PARAMETERS ###
    # NOTE: If coarse-graining, double-check your r-cut value
    # relative to your coarse-grained potentials
    # Also be sure to set the focefield parameter to None
    parameters["coarse_grain"] = [True]
    parameters["ref_distance"] = [3.3997] # Angstrom
    parameters["ref_mass"] = [15.99] # AMU
    parameters["ref_energy"] = [0.21] # kcal/mol
    # Location in polybinder.library.forcefields to look for CG specific table
    # Leave as None if the table potentials are in polybinder.library.forcefields
    parameters["cg_potentials_dir"] = [
            #"msibi-shallow",
            "msibi-deep"
    ]
    parameters["bead_mapping"] = ["ring_plus_linkage_AA"]

    ### SIMULATION PARAMETERS ###
    parameters["tau_kt"] = [0.3]
    parameters["tau_p"] = [None]
    parameters["pressure"] = [None]
    parameters["dt"] = [0.0003]
    parameters["r_cut"] = [4.0]
    parameters["e_factor"] = [1.0]
    parameters["sim_seed"] = [42]
    parameters["neighbor_list"] = ["Cell"]
    parameters["walls"] = [None]
    parameters["init_shrink_kT"] = [1.0]
    parameters["final_shrink_kT"] = [1.0]
    parameters["shrink_steps"] = [500]
    parameters["shrink_period"] = [1]
    parameters["procedure"] = ["quench"]
    parameters["gsd_write"] = [400000]
    parameters["log_write"] = [40000]

    ### Quench related parameters ###
    parameters["kT_quench"] = [
            2.0,
            #2.5
    ]
    parameters["n_steps"] = [
            #3e9,
            5e9,
    ]

    ### Anneal related parameters ###
    # List of [initial kT, final kT] Reduced Temps
    #parameters["kT_anneal"] = [
    #        [2.0, 0.5]
    #    ]
    # List of lists of number of steps
    #parameters["anneal_sequence"] = [
    #        [2e7]*10
    #]
    #parameters["schedule"] = [None]

    ### Do a few checks to save headaches later:
    if [None] not in [
			parameters["para_weight"], parameters["monomer_sequence"]
			]:
        raise ValueError(
                "You can only use one of `para_weight` and "
                "`monomer_sequence`. Set the other to [None]"
                )
    return list(parameters.keys()), list(product(*parameters.values()))

custom_job_doc = {} # add keys and values for each job document created


def main():
    project = signac.init_project("ti-ratio-relax") # Set the signac project name
    param_names, param_combinations = get_parameters()
    # Create the generate jobs
    for params in param_combinations:
        parent_statepoint = dict(zip(param_names, params))
        parent_job = project.open_job(parent_statepoint)
        parent_job.init()
        parent_job.doc.setdefault("done", False)
        parent_job.doc.setdefault("extra_runs", 0)
        try:
            parent_job.doc.setdefault(
                    "steps", int(parent_statepoint["n_steps"])
            )
        except KeyError:
            parent_job.doc.setdefault(
                    "steps", int(sum(parent_statepoint["anneal_sequence"]))
            )
            parent_job.doc.setdefault(
                    "step_sequence", parent_statepoint["anneal_sequence"]
            )
        if any(
            [parent_job.sp['Mn'], parent_job.sp['pdi'], parent_job.sp['Mw']]
        ):
            parent_job.doc.setdefault("sample_pdi", True)
        else:
            parent_job.doc.setdefault("sample_pdi", False)

    if custom_job_doc:
        for key in custom_job_doc:
            parent_job.doc.setdefault(key, custom_job_doc[key])

    project.write_statepoints()

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
