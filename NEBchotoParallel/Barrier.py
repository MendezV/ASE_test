from ase.io import read
from ase.constraints import FixAtoms
from ase.neb import NEB
from JDFTx import JDFTx
from ase.optimize import BFGS
from ase.parallel import rank, size

initial = read('initial.traj')
final = read('final.traj')


constraint = FixAtoms(mask=[atom.tag > 1 for atom in initial])
images = [initial]



for i in range(3):
    image = initial.copy()
    image.set_calculator(JDFTx(executable='mpirun -n 1 -N 1 -c 12 --exclude=node[1001-1032] /home/jfm343/JDFTXDIR/build/jdftx', pseudoDir='/home/jfm343/JDFTXDIR/build/pseudopotentials',pseudoSet='GBRV-pbe',commands={'elec-cutoff' : '20 100','kpoint-folding' : '2 2 1'}))
    image.set_constraint(constraint)
    images.append(image)
images.append(final)

neb = NEB(images, parallel=True)
neb.interpolate()
qn = BFGS(neb, trajectory='neb.traj')
qn.run(fmax=0.05)
