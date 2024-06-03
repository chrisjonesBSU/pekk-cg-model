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
import random
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
    parameters["system_type"] = ["stack"]
    parameters["molecule"] = ['PEKK']
    parameters["para_weight"] = [1.0, 0.80, 0.70, 0.60]
    parameters["monomer_sequence"] = [None]
    parameters["density"] = [0.0003]
    parameters["n_compounds"] = [[1]]
    parameters["polymer_lengths"] = [[20], [30], [10]]
    parameters["pdi"] = [None]
    parameters["Mn"] = [None]
    parameters["Mw"] = [None]
    parameters['mass_dist'] = ['weibull']
    parameters["charges"] = ["antechamber"]
    parameters["forcefield"] = ["gaff"]
    parameters["remove_hydrogens"] = [
            True,
            #False
    ]
    parameters["system_seed"] = [
            24,
            #1,
            #234,
            #84,
            #540
    ]
    parameters["box_constraints"] = [
            {"x": None, "y": None, "z": None}
	]
    parameters["kwargs"] = [
            {}
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
    # relative to your coarse-grained table potential files
    # Also be sure to set the focefield parameter to None
    parameters["coarse_grain"] = [False]
    parameters["ref_distance"] = [None] # Angstrom
    parameters["ref_mass"] = [None] # AMU
    parameters["ref_energy"] = [None] # kcal/mol
    # Location in polybinder.library.forcefields to look for CG specific table
    parameters["cg_potentials_dir"] = [None]
    parameters["cg_bead"] = [None]
    parameters["bead_mapping"] = [None]
    parameters["ekk_weight"] = [None]
    parameters["kek_weight"] = [None]
    parameters["dihedrals"] = [None]
    ### SIMULATION PARAMETERS ###
    parameters["tau_kt"] = [0.03]
    parameters["tau_p"] = [None]
    parameters["pressure"] = [None]
    parameters["couple"] = ["xyz"]
    parameters["dt"] = [0.0003]
    parameters["r_cut"] = [2.5]
    parameters["e_factor"] = [1.0]
    parameters["sim_seed"] = [42]
    parameters["neighbor_list"] = ["Tree"]
    parameters["walls"] = [None]
    parameters["init_shrink_kT"] = [None]
    parameters["final_shrink_kT"] = [None]
    parameters["shrink_steps"] = [0]
    parameters["shrink_period"] = [None]
    parameters["procedure"] = [
            "quench",
            #"anneal"
        ]
    parameters["num_gsd_frames"] = [2000]
    parameters["num_log_lines"] = [10000]

    ### Quench related parameters ###
    parameters["kT_quench"] = [
            4.0,
            4.5,
            5.0,
            5.5,
            6.0,
            6.5,
            7.0,
            7.5,
            8.0,
            8.5
    ]
    parameters["n_steps"] = [5e7]
    return list(parameters.keys()), list(product(*parameters.values()))


def main():
    project = signac.init_project("rg-re-pl") # Set the signac project name
    param_names, param_combinations = get_parameters()
    # Create the generate jobs
    for params in param_combinations:
        parent_statepoint = dict(zip(param_names, params))
        parent_job = project.open_job(parent_statepoint)
        parent_job.init()
        parent_job.doc.setdefault("done", False)
        parent_job.doc.setdefault(
                "steps", parent_job.sp.n_steps + parent_job.sp.shrink_steps
        )

        random.seed(parent_job.sp.system_seed)
        p_list = ["P"]*int(parent_job.sp.polymer_lengths[0]*parent_job.sp.para_weight)
        m_list = ["M"]*int(parent_job.sp.polymer_lengths[0]*(1-parent_job.sp.para_weight))
        full_list = p_list + m_list
        random.shuffle(full_list)
        parent_job.doc.setdefault("pm_sequence", "".join(full_list))

    project.write_statepoints()

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
