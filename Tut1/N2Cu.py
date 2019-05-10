from ase import Atoms
from ase.calculators.emt import EMT
from ase.constraints import FixAtoms
from ase.optimize import QuasiNewton
from ase.build import fcc111, add_adsorbate

h = 1.85
d = 1.10


from JDFTx import JDFTx
calculator = JDFTx(
    executable='srun -n 1 -N 1 -c 12 --exclude=node[1001-1032] /home/jfm343/JDFTXDIR/build/jdftx',
    pseudoDir='/home/jfm343/JDFTXDIR/build/pseudopotentials',
    pseudoSet='GBRV-pbe',
    commands={'elec-cutoff' : '20 100','kpoint-folding' : '2 2 1', 'electronic-scf':' ', 'elec-smearing': 'Fermi 0.01', 'elec-ex-corr' :'gga-PBEsol'})

slab = fcc111('Cu', size=(3, 3, 2), vacuum=10.0)

slab.set_calculator(calculator)
e_slab = slab.get_potential_energy()

molecule = Atoms('2N', positions=[(0., 0., 0.), (0., 0., d)])
molecule.set_calculator(calculator)
e_N2 = molecule.get_potential_energy()

add_adsorbate(slab, molecule, h, 'ontop')
constraint = FixAtoms(mask=[a.symbol != 'N' for a in slab])
slab.set_constraint(constraint)
dyn = QuasiNewton(slab, trajectory='N2Cu.traj')
dyn.run(fmax=0.05)

print('Adsorption energy:', e_slab + e_N2 - slab.get_potential_energy())
