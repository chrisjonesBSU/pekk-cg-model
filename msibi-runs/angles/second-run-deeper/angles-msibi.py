from msibi import MSIBI, State, Pair, Bond, Angle, Dihedral
import numpy as np
import signac
import os
import shutil


# Call the signac project that contains the single-chain, low density target UA simulation
p_path = "/home/erjank_project/chrisjones/pekk-msibi-final/learning-runs/single-chains"
project = signac.get_project(p_path)

kT = 6.5
single_chain_runs = [
    job for job in project.find_jobs(
        {"n_compounds": [1], "kT_quench": kT, "polymer_lengths": [20], "system_seed": 24}
    )
]


for idx, job in enumerate(single_chain_runs):
    n_steps = 6e6
    weight = job.sp.para_weight
    print(f"Run {idx} for para weight {weight}")

    opt = MSIBI(
        nlist="hoomd.md.nlist.tree",
        integrator="hoomd.md.integrate.nvt",
        integrator_kwargs={"tau": 0.01},
        dt=0.0003,
        gsd_period=int(n_steps/1000),
        n_steps=n_steps,
    )
    ## Create State object, and add it to the opt.states attribute
    ## Only using a single state to optimize bonded potentials
    opt.add_state(
        State(
            name="A",
            kT=job.sp.kT_quench,
            traj_file=job.fn("components.gsd"),
            alpha=0.7,
            max_frames=500
        )
    )

    ## Create Pair objects, and add them to the opt.pairs attribute
    ## For optimizing the bond-stretching potential, pair potentials are
    # "turned off" (LJ potential with epsilon=0)
    pair_path = "/home/erjank_project/chrisjones/pekk-msibi-final/msibi-runs/pairs-6.5kT/"
    pair0 = Pair(type1="E", type2="E")
    pair1 = Pair(type1="K", type2="K")
    pair2 = Pair(type1="E", type2="K")
    for pair in [pair0, pair1, pair2]:
        pair.set_from_file(
                os.path.join(
                pair_path,
                "workspace",
                "cf6f340fd41d9f9119443c10dc659a95",
                f"{pair.type1}-{pair.type2}_final.txt"
            )
        )
        opt.add_pair(pair)

    ## Create Bond objects, and add them to the opt.bonds
    bond0 = Bond(type1="E", type2="K")
    bond1 = Bond(type1="K", type2="K")
    bond0.set_harmonic(k=850, l0=1.47)
    bond1.set_harmonic(k=1450, l0=1.53)
    opt.add_bond(bond0)
    opt.add_bond(bond1)

    # Create Angle objects, and add them to opt.angles
    # Since we are optimizing angles, set quadratic pot with a guess
    angle_path = "/home/erjank_project/chrisjones/pekk-msibi-final/cg-potentials"
    angle0 = Angle(type1="E", type2="K", type3="K", head_correction_form="linear")
    angle0.set_from_file(
            os.path.join(angle_path,f"E-K-K_angle_{weight}.txt")
        )

    angle1 = Angle(type1="K", type2="E", type3="K", head_correction_form="linear")
    angle1.set_from_file(
            os.path.join(angle_path,f"K-E-K_angle_{weight}.txt")
        )

    opt.add_angle(angle0)
    opt.add_angle(angle1)

    target_path = "/home/erjank_project/chrisjones/pekk-msibi-final/learning-runs/single-chains/average_target_angle_dists/6.5kT"
    angle0_target = np.loadtxt(
            os.path.join(target_path, f"ekk_target_dist_6.5kT_{weight}_TI.txt")
    )
    angle1_target = np.loadtxt(
            os.path.join(target_path, f"kek_target_dist_6.5kT_{weight}_TI.txt")
    )
    angle0._set_target_distribution(
            state=opt.states[0], target_distribution=angle0_target
    )
    angle1._set_target_distribution(
            state=opt.states[0], target_distribution=angle1_target
    )

    dihedral0 = Dihedral(type1="E", type2="K", type3="K", type4="E")
    dihedral0.set_harmonic(k=16, d=-1, n=1, phi0=0)
    dihedral1 = Dihedral(type1="K", type2="E", type3="K", type4="K")
    dihedral1.set_harmonic(k=12, d=-1, n=1, phi0=0)
    opt.add_dihedral(dihedral0)
    opt.add_dihedral(dihedral1)

    ## Run the optimization
    opt.optimize_angles(n_iterations=8, smooth_pot=True, smoothing_window=7)

    # Set up P/M Dirs and move results
    os.mkdir(os.path.join(os.getcwd(), f"ekk_{weight}_{kT}kT"))
    for _dir in ["potentials", "states"]:
        shutil.move(
            os.path.join(os.getcwd(), _dir),
            os.path.join(os.getcwd(), f"ekk_{weight}_{kT}kT")
        )
