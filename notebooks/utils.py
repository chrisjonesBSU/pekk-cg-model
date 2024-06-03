import gsd.hoomd
import MDAnalysis as mda
from MDAnalysis.analysis import polymer
import numpy as np
from polybinderCG.coarse_grain import System
from cmeutils.sampling import equil_sample


def lj_potential(r, sigma, epsilon, m=12, n=6):
    v = 4*epsilon*((sigma/r)**m-(sigma/r)**n)
    return v


def persistence_length(gsd, start, stop):
    u = mda.Universe(gsd)
    chains = u.atoms.fragments
    backbones = [chain.select_atoms("name E K") for chain in chains]
    sorted_backbones = [polymer.sort_backbone(bb) for bb in backbones]
    pl = polymer.PersistenceLength(sorted_backbones)
    pl = pl.run(start=start, stop=stop)
    return pl


def radius_of_gyration(gsd, n_frames):
    sys = System(gsd_file=gsd, atoms_per_monomer=1)
    rg_values = []
    for i in range(0, n_frames):
        sys.update_frame(frame=-(i+1))
        rg = sys.radii_of_gyration(use_monomers=True)
        rg_values.extend(rg)
    return rg_values


def end_to_end(gsd, n_frames, mean=False):
    sys = System(gsd_file=gsd, atoms_per_monomer=1)
    re_values = []
    for i in range(0, n_frames):
        sys.update_frame(frame=-(i+1))
        re = sys.end_to_end_distances()
        re_values.extend(re)
    if mean:
        re_values = np.mean(re_values)
    return re_values


def sample_Rg(gsd_file, cut_frames):
    with gsd.hoomd.open(gsd_file) as traj:
        n_frames = len(traj)
    rg_values = radius_of_gyration(gsd_file, n_frames)
    rg_values.reverse()
    return equil_sample(np.array(rg_values[cut_frames:]))


def sample_Re(gsd_file, cut_frames):
    with gsd.hoomd.open(gsd_file) as traj:
        n_frames = len(traj)
    re_values = end_to_end(gsd_file, n_frames)
    re_values.reverse()
    return equil_sample(np.array(re_values[cut_frames:]))
