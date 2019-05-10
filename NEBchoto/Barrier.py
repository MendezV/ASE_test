from ase.io import read
from ase.constraints import FixAtoms
from JDFTx import JDFTx
from ase.neb import NEB
from ase.optimize import  BFGS, LBFGS

initial = read('initial.traj')
final = read('final.traj')

constraint = FixAtoms(mask=[atom.tag > 1 for atom in initial])

for atom in initial:
    print(atom.tag,atom.symbol)
'''
images = [initial]
for i in range(3):
    image = initial.copy()

    image.set_calculator(JDFTx(
    executable='srun -n 1 -N 1 -c 12 --exclude=node[1001-1032] /home/jfm343/JDFTXDIR/build/jdftx',
    pseudoDir='/home/jfm343/JDFTXDIR/build/pseudopotentials',
    pseudoSet='GBRV-pbe',
    commands={'elec-cutoff' : '20 100','kpoint-folding' : '2 2 1'}))
    
    image.set_constraint(constraint)
    images.append(image)

images.append(final)

neb = NEB(images)
neb.interpolate()
dyn = LBFGS(neb, trajectory='neb.traj')
dyn.run(fmax=0.05)
'''
