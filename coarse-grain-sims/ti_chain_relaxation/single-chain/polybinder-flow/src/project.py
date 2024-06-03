"""Define the project's workflow logic and operation functions.

Execute this script directly from the command line, to view your project's
status, execute operations and submit them to a cluster. See also:

    $ python src/project.py --help
"""
import signac
from flow import FlowProject, directives
from flow.environment import DefaultSlurmEnvironment
import os


class MyProject(FlowProject):
    pass


class Borah(DefaultSlurmEnvironment):
    hostname_pattern = "borah"
    template = "borah.sh"

    @classmethod
    def add_args(cls, parser):
        parser.add_argument(
            "--partition",
            default="gpu",
            help="Specify the partition to submit to."
        )


class R2(DefaultSlurmEnvironment):
    hostname_pattern = "r2"
    template = "r2.sh"

    @classmethod
    def add_args(cls, parser):
        parser.add_argument(
            "--partition",
            default="gpuq",
            help="Specify the partition to submit to."
        )


class Fry(DefaultSlurmEnvironment):
    hostname_pattern = "fry"
    template = "fry.sh"

    @classmethod
    def add_args(cls, parser):
        parser.add_argument(
            "--partition",
            default="batch",
            help="Specify the partition to submit to."
        )

# Definition of project-related labels (classification)
@MyProject.label
def sampled(job):
    return job.doc.get("done")


@MyProject.label
def initialized(job):
    if job.sp.coarse_grain == False:
        return job.isfile("init.mol2")
    else:
        return job.isfile("atomistic_gsd.gsd")


def get_gsd_file(job):
    if job.sp.signac_project and job.sp.signac_args:
        print("Restarting job from another signac workspace")
        project = signac.get_project(
            root=job.sp.signac_project, search=False
        )
        print("Found project:")
        print(project)
        if isinstance(job.sp.signac_args, signac.core.attrdict.SyncedAttrDict):
            print("-------------------------------")
            print("Restart job filter used:")
            print(job.sp.signac_args)
            print("-------------------------------")
            job_lookup = list(project.find_jobs(filter=job.sp.signac_args))
            if len(job_lookup) > 1:
                print([j.id for j in job_lookup])
                raise ValueError(
                        "The signac filter provied returned more than "
                        "1 job."
                )
            if len(job_lookup) < 1:
                raise ValueError(
                        "The signac filter provided found zero jobs."
                )
            _job = job_lookup[0]
            print("-------------------------------")
            print("Found restart job:")
            print(_job.id)
            print(_job.sp)
            print("-------------------------------")
        elif isinstance(job.sp.signac_args, str): # Find job using job ID
            print("Restart job found by job ID:")
            _job = project.open_job(id=job.sp.signac_args)
            print(f"Job ID: {_job.id}")
        restart_file = _job.fn('restart.gsd')
    elif job.sp.slab_file:
        restart_file = job.sp.restart_file
    return restart_file, _job.doc.final_timestep


@directives(executable="python -u")
@directives(ngpu=1)
@MyProject.operation
@MyProject.post(sampled)
def sample(job):
    from polybinder import simulate, system
    from polybinder.utils import base_units, unit_conversions
    import numpy as np
    import hoomd

    with job:
        print("-----------------------")
        print("JOB ID NUMBER:")
        print(job.id)
        print("-----------------------")
        print("----------------------")
        print("Creating the system...")
        print("----------------------")
        system_parms = system.System(
                density=job.sp.density,
                molecule=job.sp.molecule,
                n_compounds=list(job.sp.n_compounds),
                polymer_lengths=list(job.sp.polymer_lengths),
                para_weight=job.sp.para_weight,
                monomer_sequence=job.sp.monomer_sequence,
                sample_pdi=job.doc.sample_pdi,
                pdi=job.sp.pdi,
                Mn=job.sp.Mn,
                Mw=job.sp.Mw,
                seed=job.sp.system_seed
        )
        system = system.Initializer(
                system=system_parms,
                forcefield=job.sp.forcefield,
                charges=job.sp.charges,
                remove_hydrogens=job.sp.remove_hydrogens,
        )

        ref_values = None
        auto_scale = True

        if job.isfile("restart.gsd"): # Restarting from same workspace
            print("--------------------------------------------------")
            print("Initializing simulation from a restart.gsd file...")
            print("--------------------------------------------------")
            restart = job.fn("restart.gsd")
            init_shrink_kT = None
            final_shrink_kT = None
            shrink_steps = 0
            shrink_period = None
        else: # Not restarting job
            restart = None
            init_shrink_kT = job.sp.init_shrink_kT
            final_shrink_kT = job.sp.final_shrink_kT
            shrink_steps = job.sp.shrink_steps
            shrink_period = job.sp.shrink_period

        if job.sp.coarse_grain == True:
            print("----------------------------------------")
            print("Preparing a coarse-grained simulation...")
            print("----------------------------------------")
            system.coarse_grain_system(
                    bead_mapping=job.sp.bead_mapping,
                    use_components=True
            )
            ref_values = {
                "distance": job.sp.ref_distance,
                "energy": job.sp.ref_energy,
                "mass": job.sp.ref_mass
            }
            auto_scale = False
            cg_potentials_dir = job.sp.cg_potentials_dir
            #n_steps = job.sp.n_steps

        if job.sp.coarse_grain == False:
            system.system.save('init.mol2', overwrite=True)
            cg_potentials_dir = None

        system.stack()

        print("-------------------")
        print("System generated...")
        print("-------------------")
        print("----------------------")
        print("Starting simulation...")
        print("----------------------")
        simulation = simulate.Simulation(
                system,
                r_cut=job.sp.r_cut,
                tau_kt=job.sp.tau_kt,
                tau_p=job.sp.tau_p,
                nlist=job.sp.neighbor_list,
                wall_axis=job.sp.walls,
                dt=job.sp.dt,
                seed=job.sp.sim_seed,
                auto_scale=auto_scale,
                ref_values=ref_values,
                mode="gpu",
                gsd_write=max([int(job.doc.steps/10000), 1]),
                log_write=max([int(job.doc.steps/20000), 1]),
                restart=restart,
				cg_potentials_dir=cg_potentials_dir,
                ekk_weight=job.sp.ekk_weight,
                kek_weight=job.sp.kek_weight,
                bond_kwargs=job.sp.bonds,
                dihedral_kwargs=job.sp.dihedrals,
                nlist_exclusions=["bond", "angle"]
        )
        print("------------------------------")
        print("Simulation object generated...")
        print("------------------------------")
        hoomd.write.GSD.write(
                simulation.sim.state, filename=os.path.join(job.ws, "init.gsd")
        )
        job.doc['ref_energy'] = simulation.ref_energy
        job.doc['ref_distance'] = simulation.ref_distance
        job.doc['ref_mass'] = simulation.ref_mass
        job.doc['real_timestep'] = unit_conversions.convert_to_real_time(
				simulation.dt,
                simulation.ref_energy,
                simulation.ref_distance,
                simulation.ref_mass
        )
        job.doc['time_unit'] = 'fs'
        job.doc['steps_per_frame'] = simulation.gsd_write
        job.doc['steps_per_log'] = simulation.log_write

        print("-----------------------------")
        print("Quench simulation started...")
        print("-----------------------------")
        simulation.quench(n_steps=job.sp.n_steps, kT=job.sp.kT_quench)
        print("-----------------------------")
        print("Quench simulation finished...")
        print("-----------------------------")

        job.doc.done = True
        print("-----------------------------")
        print("Simulation finished completed")
        print("-----------------------------")

if __name__ == "__main__":
    MyProject().main()
