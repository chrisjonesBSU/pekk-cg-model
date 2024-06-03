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
def nvt_done(job):
    return job.doc.nvt_done


@MyProject.label
def sample_done(job):
    return job.doc.sample_done


@MyProject.post(nvt_done)
@MyProject.operation(
        directives={"ngpu": 1, "executable": "python -u"}, name="nvt"
)
def run_nvt(job):
    import unyt as u
    from unyt import Unit
    import cmeutils
    from cmeutils.gsd_utils import ellipsoid_gsd

    import flowermd
    from flowermd.base import Pack, Lattice, Simulation
    from flowermd.library.polymers import EllipsoidChain
    from flowermd.library.forcefields import EllipsoidForcefield
    from flowermd.utils.rigid_body import create_rigid_body
    from flowermd.utils import get_target_box_number_density

    with job:
        print("JOB ID NUMBER:")
        print(job.id)
        print("------------------------------------")
        chains = EllipsoidChain(
                lengths=job.sp.lengths,
                num_mols=job.sp.num_mols,
                lpar=job.sp.lpar,
                bead_mass=job.sp.bead_mass,
                bond_length=0.001,
        )
        #system = Lattice(
        #        molecules=chains,
        #        n=8,
        #        y=1.2,
        #        x=1.2,
        #        base_units = {
        #            "mass": job.sp.bead_mass * Unit("amu"),
        #            "length": job.sp.lpar * Unit("nm"),
        #            "energy": job.sp.epsilon * Unit("kJ/mol")
        #        }
        #)
        system = Pack(
                molecules=chains,
                density=job.sp.density * (1/u.Unit("nm**3")),
                fix_orientation=True,
                base_units = {
                    "mass": job.sp.bead_mass * Unit("amu"),
                    "length": job.sp.lpar * Unit("nm"),
                    "energy": job.sp.epsilon * Unit("kJ/mol")
                }
        )
        rigid_frame, rigid = create_rigid_body(
                system.hoomd_snapshot,
                chains.bead_constituents_types
        )

        ellipsoid_ff = EllipsoidForcefield(
                epsilon=job.sp.epsilon,
                lperp=job.sp.lperp,
                lpar=job.sp.lpar,
                r_cut=job.sp.r_cut,
                bond_k=job.sp.bond_k,
                bond_r0=job.sp.r0,
                angle_k=job.sp.theta_k,
                angle_theta0=job.sp.theta0
        )

        sim = Simulation(
                initial_state=rigid_frame,
                forcefield=ellipsoid_ff.hoomd_forces,
                rigid_constraint=rigid,
                dt=job.sp.dt,
                seed=job.sp.seed,
                reference_values=system.reference_values,
                gsd_write_freq=job.sp.gsd_write_freq,
                log_write_freq=job.sp.log_write_freq,
                gsd_file_name=job.fn("trajectory.gsd"),
                log_file_name=job.fn("data.txt"),
        )
        sim.pickle_forcefield(job.fn("forcefield.pickle"))
        sim.save_restart_gsd(job.fn("init.gsd"))
        # Store unit information in job doc
        tau_kT = sim.dt * job.sp.tau_kT
        job.doc.tau_kT = tau_kT
        job.doc.ref_mass = sim.reference_mass.to("amu").value
        job.doc.ref_mass_units = "amu"
        job.doc.ref_energy = sim.reference_energy.to("kJ/mol").value
        job.doc.ref_energy_units = "kJ/mol"
        job.doc.ref_length = sim.reference_length.to("nm").value
        job.doc.ref_length_units = "nm"
        job.doc.real_time_step = sim.real_timestep.to("fs").value
        job.doc.real_time_units = "fs"
        n_beads = 0
        for n, l in zip(job.sp.lengths, job.sp.num_mols):
            n_beads += (n * l)
        job.doc.n_beads = int(n_beads)
        # Set up stuff for shrinking volume step
        print("Running shrink step.")
        shrink_kT_ramp = sim.temperature_ramp(
                n_steps=job.sp.shrink_n_steps,
                kT_start=job.sp.shrink_kT,
                kT_final=job.sp.kT
        )

        sigma = 1 * Unit("nm")
        density = job.sp.density / (sigma**3)
        target_box = get_target_box_number_density(
                density=density,
                n_beads=job.doc.n_beads
        )
        print(target_box)
        job.doc.target_box = target_box
        #sim.run_update_volume(
        #        final_box_lengths=target_box,
        #        n_steps=job.sp.shrink_n_steps,
        #        period=job.sp.shrink_period,
        #        tau_kt=tau_kT,
        #        kT=shrink_kT_ramp
        #)
        print("Shrink step finished.")
        print("Running simulation.")
        sim.run_NVT(kT=job.sp.kT, n_steps=job.sp.n_steps, tau_kt=tau_kT)
        sim.save_restart_gsd(job.fn("restart.gsd"))
        sim.flush_writers()
        job.doc.nvt_done = True
        print("Simulation finished.")
        ellipsoid_gsd(
                job.fn("trajectory.gsd"),
                job.fn("ellipsoid-trajectory.gsd"),
                lpar=job.sp.lpar,
                lperp=job.sp.lperp,
        )


@MyProject.pre(nvt_done)
@MyProject.post(sample_done)
@MyProject.operation(
        directives={"ngpu": 0, "executable": "python -u"}, name="sample"
)
def sample(job):
    # Add package imports here
    with job:
        print("JOB ID NUMBER:")
        print(job.id)
        print("------------------------------------")
        # Add your script here
        job.doc.sample_done = True


if __name__ == "__main__":
    MyProject(environment=Fry).main()
