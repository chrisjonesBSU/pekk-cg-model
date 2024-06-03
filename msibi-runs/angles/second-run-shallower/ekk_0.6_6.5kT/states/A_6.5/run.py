
import hoomd
import hoomd.md
from hoomd.init import read_gsd

hoomd.context.initialize("")
system = read_gsd("/home/erjank_project/chrisjones/pekk-msibi-final/learning-runs/single-chains/workspace/0f2d30692a6cdb8fd90ae3f6aec03a03/components.gsd", frame=-1, time_step=0)

nl = hoomd.md.nlist.tree()
nl.reset_exclusions(exclusions=['1-2', '1-3'])

table=hoomd.md.pair.table(width=101,nlist=nl)
table.set_from_file('E', 'E', filename='/home/erjank_project/chrisjones/pekk-msibi-final/msibi-runs/pairs-6.5kT/workspace/b0a23d31f30e96810cc01f23a4fe909f/E-E_final.txt')
table.set_from_file('K', 'K', filename='/home/erjank_project/chrisjones/pekk-msibi-final/msibi-runs/pairs-6.5kT/workspace/b0a23d31f30e96810cc01f23a4fe909f/K-K_final.txt')
table.set_from_file('E', 'K', filename='/home/erjank_project/chrisjones/pekk-msibi-final/msibi-runs/pairs-6.5kT/workspace/b0a23d31f30e96810cc01f23a4fe909f/E-K_final.txt')
harmonic_bond = hoomd.md.bond.harmonic()
harmonic_bond.bond_coeff.set('E-K', k=850, r0=1.47)
harmonic_bond.bond_coeff.set('K-K', k=1450, r0=1.53)
atable = hoomd.md.angle.table(width=200)
atable.set_from_file('E-K-K', '/home/erjank_project/chrisjones/pekk-msibi-final/msibi-runs/angles/second-run-shallower/potentials/angle_pot.E-K-K')
atable.set_from_file('K-E-K', '/home/erjank_project/chrisjones/pekk-msibi-final/msibi-runs/angles/second-run-shallower/potentials/angle_pot.K-E-K')
harmonic_dihedral = hoomd.md.dihedral.harmonic()
harmonic_dihedral.dihedral_coeff.set('E-K-K-E', k=16, d=-1, n=1, phi0=0)
harmonic_dihedral.dihedral_coeff.set('K-E-K-K', k=12, d=-1, n=1, phi0=0)

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

