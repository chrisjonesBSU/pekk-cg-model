
import hoomd
import hoomd.md
from hoomd.init import read_gsd

hoomd.context.initialize("")
system = read_gsd("/home/erjank_project/chrisjones/pekk-msibi-final/learning-runs/single-chains/workspace/b83dce273531102fa2c624e593dddb4c/components.gsd", frame=-1, time_step=0)

nl = hoomd.md.nlist.tree()
nl.reset_exclusions(exclusions=['1-2', '1-3'])

lj = hoomd.md.pair.lj(nlist=nl, r_cut=0)
lj.pair_coeff.set('E', 'E', epsilon=0, sigma=1, r_cut=0)
lj.pair_coeff.set('K', 'K', epsilon=0, sigma=1, r_cut=0)
lj.pair_coeff.set('E', 'K', epsilon=0, sigma=1, r_cut=0)
harmonic_bond = hoomd.md.bond.harmonic()
harmonic_bond.bond_coeff.set('E-K', k=850, r0=1.47)
harmonic_bond.bond_coeff.set('K-K', k=1450, r0=1.53)
atable = hoomd.md.angle.table(width=200)
atable.set_from_file('E-K-K', '/home/erjank_project/chrisjones/pekk-msibi-final/msibi-runs/angles-harmonic-bonds-6.5kT/potentials/angle_pot.E-K-K.txt')
atable.set_from_file('K-E-K', '/home/erjank_project/chrisjones/pekk-msibi-final/msibi-runs/angles-harmonic-bonds-6.5kT/potentials/angle_pot.K-E-K.txt')

_all = hoomd.group.all()
hoomd.md.integrate.mode_standard(0.0003)
integrator_kwargs = {'tau': 0.01, 'kT': 6.5}
integrator = hoomd.md.integrate.nvt(group=_all, **integrator_kwargs)


hoomd.dump.gsd(
    filename="query.gsd",
    group=_all,
    period=6000,
    overwrite=True,
    dynamic=["momentum"]
    )

hoomd.run(6000000)

