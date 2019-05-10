from JDFTx import JDFTx
from ase.io import read
from ase.constraints import FixAtoms
from ase.neb import NEB
from ase.optimize import BFGS,LBFGS

initial = read('initial.traj')
final = read('final.traj')

c2 = FixAtoms(indices=[4,5,6])
images = [initial]
for i in range(3):
    image = initial.copy()
    image.set_calculator(JDFTx(executable='srun -n 1 -N 1 -c 12 --exclude=node[1001-1032]  /home/jfm343/JDFTXDIR/build/jdftx',pseudoDir='/home/jfm343/JDFTXDIR/build/pseudopotentials',pseudoSet='GBRV-pbe',commands={'elec-cutoff' : '20 100','kpoint': '.5 0 0','kpoint-folding' : '2 2 2', 'electronic-scf':' ', 'elec-smearing': 'Fermi 0.01'}))
    images.append(image)
    image.set_constraint(c2) ## so that beads are able to move in intermediate steps

images.append(final)

neb = NEB(images, parallel=True)
neb.interpolate()
qn = LBFGS(neb, trajectory='neb.traj')
qn.run(fmax=0.05)

