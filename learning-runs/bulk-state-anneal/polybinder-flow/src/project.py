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
            default="shortgpu",
            help="Specify the partition to submit to."
        )


class R2(DefaultSlurmEnvironment):
    hostname_pattern = "r2"
    template = "r2.sh"

    @classmethod
    def add_args(cls, parser):
        parser.add_argument(
            "--partition",
            default="shortgpuq",
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
    from polybinder.utils import unit_conversions
    from polybinder.utils.base_units import base_units
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

        if job.isfile("pmd_structures/parmed.pickle"):
            parmed_dir = os.path.join(job.ws, "pmd_structures")
            print("----------------------")
            print("Found an existing Parmed pickle file...")
            print("----------------------")
        else:
            parmed_dir = "pmd_structures"

        # Set up system parameters
        if job.sp.system_type != "interface":
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

            job.doc['num_para'] = system_parms.para
            job.doc['num_meta'] = system_parms.meta
            job.doc['num_compounds'] = system_parms.n_compounds
            job.doc['polymer_lengths'] = system_parms.polymer_lengths
            job.doc["chain_sequences"] = system_parms.molecule_sequences

            # Call the Initializer class
            system = system.Initializer(
                    system=system_parms,
                    forcefield=job.sp.forcefield,
                    charges=job.sp.charges,
                    remove_hydrogens=job.sp.remove_hydrogens,
                    parmed_dir=parmed_dir
            )

            job.doc["total_mass_amu"] = system.system_mass
            job.doc["total_mass_g"] = system.system_mass * base_units()["amu_to_kg"]

            # Coarse-grain the system if needed
            if job.sp.coarse_grain:
                print("----------------------------------------")
                print("Preparing a coarse-grained simulation...")
                print("----------------------------------------")
                if job.sp.cg_bead == "components":
                    use_monomers=False
                    use_components=True
                elif job.sp.cg_bead == "monomers":
                    use_monomers = True
                    use_components = False
                system.coarse_grain_system(
                        use_monomers=use_monomers,
                        use_components=use_components,
                        bead_mapping=job.sp.bead_mapping
                )
                ref_values = {
                        "distance": job.sp.ref_distance,
                        "energy": job.sp.ref_energy,
                        "mass": job.sp.ref_mass
                }
                auto_scale = False
                cg_potentials_dir = job.sp.cg_potentials_dir
            # Not using a coarse-grain system
            else:
                ref_values = None
                auto_scale = True
                cg_potentials_dir = None

            # Call the correct system builder function
            print("----------------------------------------")
            print("Building the system...")
            print("----------------------------------------")
            if job.sp.system_type == "pack":
                system.pack(**job.sp.kwargs)
            elif job.sp.system_type == "stack":
                system.stack(**job.sp.kwargs)
            elif job.sp.system_type == "crystal":
                system.crystal(**job.sp.kwargs)

            job.doc["n_particles"] = len(system.system.atoms)

            # Override the default target box
            if any(list(job.sp.box_constraints.values())):
                system.set_target_box(
                        job.sp.box_constraints["x"],
                        job.sp.box_constraints["y"],
                        job.sp.box_constraints["z"]
                )

            job.doc["target_box_nm"] = system.target_box
            job.doc["target_volume_nm"] = np.prod(system.target_box)
            job.doc["target_volume_cm"] = (
                    job.doc.target_volume_nm / base_units()["cm_to_nm"]**3
            )

            # Restarting job from within the same workspace
            if job.isfile("restart.gsd"):
                print("--------------------------------------------------")
                print("Initializing simulation from a restart.gsd file...")
                print("--------------------------------------------------")
                restart = job.fn("restart.gsd")
                init_shrink_kT = None
                final_shrink_kT = None
                shrink_steps = 0
                shrink_period = None
            elif any(
                    [
                        all([job.sp.signac_project, job.sp.signac_args]),
                        job.sp.restart_file
                    ]
            ):
                print("--------------------------------------------------")
                print("Initializing simulation from a restart.gsd file...")
                print("--------------------------------------------------")
                restart, last_n_steps = get_gsd_file(job)
                print(f"Initializing from {restart}")
                init_shrink_kT = None
                final_shrink_kT = None
                shrink_steps = 0
                shrink_period = None
            # Not initializing from a restart.gsd file
            else:
                restart = None
                init_shrink_kT = job.sp.init_shrink_kT
                final_shrink_kT = job.sp.final_shrink_kT
                shrink_steps = job.sp.shrink_steps
                shrink_period = job.sp.shrink_period

        elif job.sp.system_type == "interface":
            print("-------------------------")
            print("Creating the interface...")
            print("-------------------------")
            slab_files = []
            ref_distances = []
            if job.doc.use_signac:
                signac_args = []
                if isinstance(job.sp.signac_args, list):
                    slab_1_arg = job.sp.signac_args[0]
                    signac_args.append(slab_1_arg)
                    if len(job.sp.signac_args) == 2:
                        slab_2_arg = job.sp.signac_args[1]
                        signac_args.append(slab_2_args)
                elif not isinstance(job.sp.signac_args, list):
                    signac_args.append(job.sp.signac_args)

                project = signac.get_project(
                        root=job.sp.signac_project, search=True
                )
                for arg in signac_args:
                    if isinstance(arg, dict):
                        _job = list(project.find_jobs(filter=arg))[0]
                        slab_files.append(_job.fn('restart.gsd'))
                        ref_distances.append(_job.doc['ref_distance']/10)
                    elif isinstance(arg, str): # Find job using job ID
                        _job = project.open_job(id=arg)
                        slab_files.append(_job.fn('restart.gsd'))
                        ref_distances.append(_job.doc['ref_distance']/10)
            else:
                slab_files.append(job.sp.slab_file)
                ref_distances.append(job.sp.reference_distance)

            system = system.Interface(
					slabs=slab_files,
                    ref_distance=ref_distances[0],
                    gap=job.sp.interface_gap,
					weld_axis=job.sp.weld_axis,
            )

            job.doc['slab_ref_distances'] = system.ref_distance
            init_shrink_kT = None
            final_shrink_kT = None
            shrink_steps = 0
            shrink_period = None

        print("----------------------")
        print("System generated...")
        print("----------------------")
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
                gsd_write=max([int(job.doc.steps/job.sp.num_gsd_frames), 1]),
                log_write=max([int(job.doc.steps/job.sp.num_log_lines), 1]),
                restart=restart,
				cg_potentials_dir=cg_potentials_dir,
                ekk_weight=job.sp.ekk_weight,
                kek_weight=job.sp.kek_weight,
                dihedral_kwargs=job.sp.dihedrals,
        )
        target_box = np.copy(simulation.target_box)
        simulation.target_box = target_box * 1.25
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
        print("----------------------------")
        print("Running shrink simulation...")
        print("----------------------------")
        simulation.shrink(
                kT_init=init_shrink_kT,
                kT_final=job.sp.kT_quench,
                period=shrink_period,
                n_steps=shrink_steps,
                tree_nlist=False
        )
        print("-----------------------------")
        print("Shrink simulation finished...")
        print("-----------------------------")
        print("Quench at low density...")
        simulation.quench(
                kT=job.sp.kT_quench,
                n_steps=2e8,
                pressure=None
        )
        print("Quench finished...")
        simulation.target_box = target_box * 0.90
        print("Compressing...")
        simulation.shrink(
                kT_init=job.sp.kT_quench,
                kT_final=job.sp.kT_quench,
                period=int(1e3),
                n_steps=1e7,
                tree_nlist=False
        )
        print("Compressing finished...")
        simulation.target_box = target_box
        print("Expanding...")
        simulation.shrink(
                kT_init=job.sp.kT_quench,
                kT_final=job.sp.kT_quench,
                period=int(1e3),
                n_steps=5e7,
                tree_nlist=False
        )
        print("Expanding finished...")
        # Run a quench simulation (NVT or NPT)
        if job.sp.procedure == "quench":
            job.doc['T_SI'] = unit_conversions.kelvin_from_reduced(
                job.sp.kT_quench, simulation.ref_energy
            )
            job.doc['T_unit'] = 'K'
            print("----------------------------")
            print("Running quench simulation...")
            print("----------------------------")
            simulation.quench(
                n_steps=job.sp.n_steps,
                kT=job.sp.kT_quench,
                pressure=job.sp.pressure
            )
            print("-----------------------------")
            print("Quench simulation finished...")
            print("-----------------------------")

        elif job.sp.procedure == "anneal":
            print("----------------------------")
            print("Running anneal simulation...")
            print("----------------------------")
            step_sequence = job.sp.anneal_sequence
            if not job.sp.schedule:
                kT_list = np.linspace(
                        job.sp.kT_anneal[0],
                        job.sp.kT_anneal[1],
                        len(job.sp.anneal_sequence),
                )
                kT_SI = [
                        unit_conversions.kelvin_from_reduced(
                            kT, simulation.ref_energy
                        ) for kT in kT_list
                ]
                job.doc['T_SI'] = kT_SI
                job.doc['T_unit'] = 'K'

            simulation.anneal(
                kT_init=job.sp.kT_anneal[0],
                kT_final=job.sp.kT_anneal[1],
                pressure=job.sp.pressure,
                step_sequence=step_sequence,
                schedule=job.sp.schedule,
            )
            print("-----------------------------")
            print("Anneal simulation finished...")
            print("-----------------------------")

        job.doc["final_timestep"] = simulation.sim.timestep
        job.doc["simulation_duration_ns"] = job.doc.final_timestep * 1e-6
        job.doc["done"] = True

        print("-----------------------------")
        print("Simulation finished completed")
        print("-----------------------------")

if __name__ == "__main__":
    MyProject().main()
