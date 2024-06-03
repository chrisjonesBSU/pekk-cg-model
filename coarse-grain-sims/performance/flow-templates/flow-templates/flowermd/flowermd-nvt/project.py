"""Define the project's workflow logic and operation functions.

Execute this script directly from the command line, to view your project's
status, execute operations and submit them to a cluster. See also:

    $ python src/project.py --help
"""
import signac
from flow import FlowProject, directives
from flow.environment import DefaultSlurmEnvironment
import os
import shutil


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
def seed_done(job):
    seed_job = [j for j in job.project.find_jobs(
            {
                "n_duplicates": None,
                "coarse_grain": job.sp.coarse_grain,
                "remove_hydrogens": job.sp.remove_hydrogens 
            }
    )][0]
    return seed_job.isfile("seed-restart.gsd")


@MyProject.label
def sim_done(job):
    return job.doc.sim_done


@MyProject.label
def seed_system(job):
    return job.sp.n_duplicates is None


@MyProject.label
def not_seed_system(job):
    return job.sp.n_duplicates is not None


def duplicate_seed(job):
    import gsd
    from flowermd.modules.welding import Interface
    seed_job = [j for j in job.project.find_jobs(
            {
                "n_duplicates": None,
                "coarse_grain": job.sp.coarse_grain,
                "remove_hydrogens": job.sp.remove_hydrogens 
            }
    )][0]
    shutil.copy(
            seed_job.fn("seed-restart.gsd"),
            job.fn("seed-restart.gsd")
    )
    shutil.copy(
            seed_job.fn("forcefield.pickle"),
            job.fn("forcefield.pickle")
    )
    n_dups = 0
    while n_dups != job.sp.n_duplicates:
        for i, axis in enumerate([(1,0,0), (0,1,0), (0,0,1)]):
            if n_dups == job.sp.n_duplicates:
                break
            dup = Interface(
                    gsd_files=[job.fn("seed-restart.gsd")],
                    interface_axis=axis,
                    gap=8,
                    wall_sigma=3,
                    remove_void_particles=False
            )
            with gsd.hoomd.open(job.fn("seed-restart.gsd"), "w") as traj:
                traj.append(dup.hoomd_snapshot)
            n_dups += 1

    with gsd.hoomd.open(job.fn("seed-restart.gsd"), "r") as traj:
        snap = traj[0]
        snap.configuration.box[0:3] *= 1.2
        with gsd.hoomd.open(job.fn("init.gsd"), "w") as new_traj:
            new_traj.append(snap)
        job.doc.N = traj[0].particles.N

    job.doc.tau_kT = seed_job.doc.tau_kT
    job.doc.ref_mass = seed_job.doc.ref_mass 
    job.doc.ref_mass_units = "amu"
    job.doc.ref_energy = seed_job.doc.ref_mass 
    job.doc.ref_energy_units = "kJ/mol"
    job.doc.ref_length = seed_job.doc.ref_length 
    job.doc.ref_length_units = "nm"


@MyProject.pre(seed_system)
@MyProject.post(seed_done)
@MyProject.operation(
        directives={"ngpu": 1, "executable": "python -u"}, name="seed"
)
def make_seed(job):
    import unyt as u
    from unyt import Unit
    import flowermd
    from flowermd.base.system import Pack
    from flowermd.base.simulation import Simulation
    from flowermd.library import PEKK_para, GAFF, PEKK_CG_FF 
    from flowermd.utils import get_target_box_mass_density

    with job:
        print("JOB ID NUMBER:")
        print(job.id)
        print("------------------------------------")

        gsd_path = job.fn("seed.gsd")
        log_path = job.fn("seed-log.txt")
    
        molecules = PEKK_para(lengths=job.sp.lengths, num_mols=job.sp.num_mols)
        if job.sp.coarse_grain:
            print("Coarse-graining the system")
            molecules.coarse_grain(
                        beads={
                            "E": "c1cc(O)ccc1",
                            "K": "c1ccc(C=O)cc1"
                        }
            )
            system = Pack(molecules=molecules, density=job.sp.density) 
            system.reference_length = 0.3399669 * Unit("nm")
            system.reference_mass = 15.99 * Unit("amu")
            system.reference_energy = 1 * Unit("kJ/mol")

            ff = PEKK_CG_FF(TI_ratio=1.0)
            sim = Simulation(
                    initial_state=system.hoomd_snapshot,
                    forcefield = ff.hoomd_forces,
                    reference_values=system.reference_values,
                    dt=job.sp.dt,
                    gsd_write_freq=job.sp.gsd_write_freq,
                    gsd_file_name=gsd_path,
                    log_write_freq=job.sp.log_write_freq,
                    log_file_name=log_path,
            )
        else:
            system = Pack(molecules=molecules, density=job.sp.density) 
            system.apply_forcefield(
                    force_field=GAFF(),
                    r_cut=job.sp.r_cut,
                    auto_scale=job.sp.auto_scale,
                    scale_charges=True,
                    remove_charges=job.sp.remove_charges,
                    remove_hydrogens=job.sp.remove_hydrogens,
                    pppm_resolution=job.sp.pppm_resolution,
                    pppm_order=job.sp.pppm_order
            )
            sim = Simulation.from_system(
                    system=system,
                    dt=job.sp.dt,
                    gsd_write_freq=job.sp.gsd_write_freq,
                    gsd_file_name=gsd_path,
                    log_write_freq=job.sp.log_write_freq,
                    log_file_name=log_path,
            )

        job.doc.ref_mass = sim.reference_mass.to("amu").value
        job.doc.ref_mass_units = "amu"
        job.doc.ref_energy = sim.reference_energy.to("kJ/mol").value
        job.doc.ref_energy_units = "kJ/mol"
        job.doc.ref_length = sim.reference_length.to("nm").value
        job.doc.ref_length_units = "nm"

        

        sim.add_walls(wall_axis=(1, 0, 0), sigma=0.5, epsilon=1.0, r_cut=1.12)
        sim.add_walls(wall_axis=(0, 1, 0), sigma=0.5, epsilon=1.0, r_cut=1.12)
        sim.add_walls(wall_axis=(0, 0, 1), sigma=0.5, epsilon=1.0, r_cut=1.12)

        sim.pickle_forcefield(job.fn("forcefield.pickle"))
        # Store unit information in job doc
        tau_kT = sim.dt * job.sp.tau_kT
        job.doc.tau_kT = tau_kT
        # Set up stuff for shrinking volume step
        print("Running shrink step.")
        shrink_kT_ramp = sim.temperature_ramp(
                n_steps=job.sp.shrink_n_steps,
                kT_start=job.sp.shrink_kT,
                kT_final=job.sp.kT
        )
        target_box = get_target_box_mass_density(
                mass=system.mass.to("g"),
                density=job.sp.density * (Unit("g/cm**3") / 4)
        )
        print("TARGET BOX")
        print(target_box)
        sim.run_update_volume(
                final_box_lengths=target_box,
                n_steps=job.sp.shrink_n_steps,
                period=job.sp.shrink_period,
                tau_kt=tau_kT,
                kT=shrink_kT_ramp
        )
        print("Shrink step finished.")
        print("Running simulation.")
        sim.run_NVT(kT=9.0, n_steps=1e6, tau_kt=tau_kT)
        sim.save_restart_gsd(job.fn("seed-restart.gsd"))
        print("Seed simulation finished.")


@MyProject.pre(not_seed_system)
@MyProject.pre(seed_done)
@MyProject.post(sim_done)
@MyProject.operation(
        directives={"ngpu": 1, "executable": "python -u"}, name="nvt"
)
def run_big_system(job):
    import hoomd
    import unyt as u
    from unyt import Unit
    import flowermd
    from flowermd.base.system import Pack
    from flowermd.base.simulation import Simulation
    from flowermd.library import PEKK_para, GAFF 
    from flowermd.utils import get_target_box_mass_density
    import pickle

    with job:
        print("JOB ID NUMBER:")
        print(job.id)
        print("------------------------------------")
        duplicate_seed(job)
        forces = []
        with open(job.fn("forcefield.pickle"), "rb") as f:
            hoomd_ff = pickle.load(f)
            print("------------------------------------")
            print(hoomd_ff)
            for f in hoomd_ff:
                if isinstance(f, hoomd.md.external.wall.LJ):
                    pass
                else:
                    forces.append(f)
            print(forces)
            print("------------------------------------")

        refs_dict = {
            "length": job.doc.ref_length * Unit(job.doc.ref_length_units),
            "energy": job.doc.ref_energy * Unit(job.doc.ref_energy_units),
            "mass": job.doc.ref_mass * Unit(job.doc.ref_mass_units),
        } 
        
        sim = Simulation(
                initial_state=job.fn("init.gsd"),
                forcefield=forces,
                reference_values=refs_dict,
        )
        target_box = get_target_box_mass_density(
                mass=sim.mass.to("g"),
                density=job.sp.density * Unit("g/cm**3")
        )
        sim.run_update_volume(
                final_box_lengths=target_box,
                n_steps=1e5,
                period=1,
                tau_kt=job.doc.tau_kT,
                kT=8.0
        )
        print("Shrink step finished.")
        print("Running simulation.")
        sim.run_NVT(kT=6.0, n_steps=1e5, tau_kt=job.doc.tau_kT)
        job.doc.sim_done = True


if __name__ == "__main__":
    MyProject(environment=Fry).main()
