from ase.build import fcc100, add_adsorbate
from ase.constraints import FixAtoms
from JDFTx import JDFTx
from ase.optimize import BFGS, LBFGS

# 2x2-Al(001) surface with 3 layers and an
# Au atom adsorbed in a hollow site:
slab = fcc100('Al', size=(2, 2, 3))
add_adsorbate(slab, 'Au', 1.7, 'hollow')
slab.center(axis=2, vacuum=4.0)

# Make sure the structure is correct:
#view(slab)

# Fix second and third layers:
mask = [atom.tag > 1 for atom in slab]
#print(mask)
slab.set_constraint(FixAtoms(mask=mask))

# Use JDFTx potential:
'''
calculator = JDFTx(
    executable='/home/jfm343/JDFTXDIR/build/jdftx',
    pseudoDir='/home/jfm343/JDFTXDIR/build/pseudopotentials',
    pseudoSet='GBRV-pbe',
    commands={'elec-cutoff' : '20 100','kpoint-folding' : '6 6 1'})
'''
calculator = JDFTx(
    executable='srun -n 1 -N 1 -c 12 --exclude=node[1001-1032] /home/jfm343/JDFTXDIR/build/jdftx',
    pseudoDir='/home/jfm343/JDFTXDIR/build/pseudopotentials',
    pseudoSet='GBRV-pbe',
    commands={'elec-cutoff' : '20 100','kpoint-folding' : '2 2 1'})
slab.set_calculator(calculator)

# Initial state:
dyn = LBFGS(slab, trajectory='initial.traj')
dyn.run(fmax=0.05)


# Final state:
slab[-1].x += slab.get_cell()[0, 0] / 2
dyn = LBFGS(slab, trajectory='final.traj')
dyn.run(fmax=0.05)

calculator.clean()
