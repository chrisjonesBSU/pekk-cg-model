"""Define the project's workflow logic and operation functions.

Execute this script directly from the command line, to view your project's
status, execute operations and submit them to a cluster. See also:

    $ python src/project.py --help
"""
import signac
from flow import FlowProject, directives
from flow.environment import DefaultSlurmEnvironment
import os
import unyt
from unyt import Unit


class PEKK_Weld(FlowProject):
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
            default="v100," "batch",
            help="Specify the partition to submit to."
        )


#@MyProject.label
@PEKK_Weld.label
def initial_run_done(job):
    return job.doc.runs >= 1


#@MyProject.label
@PEKK_Weld.label
def equilibrated(job):
    return job.doc.equilibrated


def get_ref_values(job):
    ref_length = job.doc.ref_length * Unit(job.doc.ref_length_unit)
    ref_mass = job.doc.ref_mass * Unit(job.doc.ref_mass_unit)
    ref_energy = job.doc.ref_energy * Unit(job.doc.ref_energy_unit)
    ref_values_dict = {
            "length": ref_length,
            "mass": ref_mass,
            "energy": ref_energy
    }
    return ref_values_dict

def cg_mapping():
    return {"E": "c1ccc(O)cc1" ,"K": "c1ccc(C=O)cc1"}


@PEKK_Weld.post(initial_run_done)
@PEKK_Weld.operation(
        directives={"ngpu": 1, "executable": "python -u"}, name="initial-run"
)
def initial_run(job):
    """Run a bulk slab simulation; equilibrate in NVT"""
    from flowermd.base import Pack, Simulation
    from flowermd.library import PEKK_para, PEKK_CG_FF
    from flowermd.utils import get_target_box_mass_density
    import unyt as u
    from unyt import Unit

    with job:
        print("------------------------------------")
        print("JOB ID NUMBER:")
        print(job.id)
        print("------------------------------------")
        # Store reference units and values
        job.doc.ref_mass = 15.99
        job.doc.ref_mass_unit = "amu"
        job.doc.ref_energy = 0.87864
        job.doc.ref_energy_unit = "kJ/mol"
        job.doc.ref_length = 0.33996
        job.doc.ref_length_unit = "nm"

        mapping = cg_mapping()
        pekk = PEKK_para(num_mols=job.sp.num_mols, lengths=job.sp.lengths)
        pekk.coarse_grain(beads=mapping)
        system = Pack(
                molecules=pekk,
                density=job.sp.density,
                base_units=get_ref_values(job),
                edge=1.0,
                overlap=1.0,
                packing_expand_factor=7
        )
        job.doc.total_mass_amu = system.mass
        pekk_cg_ff = PEKK_CG_FF(TI_ratio=job.sp.ti_ratio)
        gsd_path = job.fn(f"trajectory{job.doc.runs + 1}.gsd")
        log_path = job.fn(f"log{job.doc.runs + 1}.txt")
        hoomd_ff = pekk_cg_ff.hoomd_forces
        sim = Simulation(
                initial_state=system.hoomd_snapshot,
                forcefield=hoomd_ff,
                reference_values=get_ref_values(job),
                dt=job.sp.dt,
                gsd_write_freq=job.sp.gsd_write_freq,
                gsd_file_name=gsd_path,
                log_write_freq=job.sp.log_write_freq,
                log_file_name=log_path,
                seed=job.sp.sim_seed
        )
        sim.pickle_forcefield(job.fn("forcefield.pickle"))
        # Store more unit information in job doc
        target_box = get_target_box_mass_density(
                mass=system.mass.to("g"),
                density=job.sp.density * Unit("g/cm**3")
        )
        tau_kT = job.sp.dt * job.sp.tau_kT
        job.doc.tau_kT = tau_kT
        job.doc.real_time_step = sim.real_timestep.to("fs").value
        job.doc.real_time_units = "fs"

        # Set up stuff for shrinking volume step
        print("Running shrink step.")
        shrink_kT_ramp = sim.temperature_ramp(
                n_steps=job.sp.shrink_n_steps,
                kT_start=job.sp.shrink_kT,
                kT_final=job.sp.kT
        )
        # Shrink
        sim.run_update_volume(
                final_box_lengths=target_box,
                n_steps=job.sp.shrink_n_steps,
                period=int(job.sp.shrink_period),
                tau_kt=job.doc.tau_kT,
                kT=shrink_kT_ramp
        )
        print("Shrinking finished.")
        print("Running NVT simulation.")
        sim.run_NVT(
                n_steps=job.sp.n_steps, kT=job.sp.kT, tau_kt=job.doc.tau_kT
        )
        sim.save_restart_gsd(job.fn("restart.gsd"))
        job.doc.runs += 1
        print("Simulation finished.")


@PEKK_Weld.pre(initial_run_done)
@PEKK_Weld.post(equilibrated)
@PEKK_Weld.operation(
        directives={"ngpu": 1, "executable": "python -u"},
        name="run-longer"
)
def run_nvt_longer(job):
    from flowermd.base import Pack, Simulation
    from flowermd.library import PEKK_para, PEKK_CG_FF
    from flowermd.utils import get_target_box_mass_density
    import unyt as u
    from unyt import Unit
    import pickle
    with job:
        print("------------------------------------")
        print("JOB ID NUMBER:")
        print(job.id)
        print("------------------------------------")
        print("Restarting NVT simulation...")
        with open(job.fn("forcefield.pickle"), "rb") as f:
            ff = pickle.load(f)

        gsd_path = job.fn(f"trajectory{job.doc.runs + 1}.gsd")
        log_path = job.fn(f"log{job.doc.runs + 1}.txt")
        ref_values = get_ref_values(job)
        sim = Simulation(
                initial_state=job.fn("restart.gsd"),
                forcefield=ff,
                reference_values=ref_values,
                dt=job.sp.dt,
                gsd_write_freq=job.sp.gsd_write_freq,
                gsd_file_name=gsd_path,
                log_write_freq=job.sp.log_write_freq,
                log_file_name=log_path,
                seed=job.sp.sim_seed,
        )
        sim.run_NVT(n_steps=1e8, kT=job.sp.kT, tau_kt=job.doc.tau_kT)
        sim.save_restart_gsd(job.fn("restart.gsd"))
        job.doc.runs += 1
        print("Simulation finished.")

if __name__ == "__main__":
    PEKK_Weld(environment=Fry).main()
