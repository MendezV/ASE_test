from ase.io import read
from ase.constraints import FixAtoms
from ase.neb import NEB
from ase.optimize import BFGS
from JDFTx import JDFTx

initial = read('initial.traj')
final = read('final.traj')

constraint = FixAtoms(mask=[0, 0, 0, 0, 1, 1, 1])

calculator = JDFTx(
    executable='srun -n 1 -N 1 -c 12 --exclude=node[1001-1032] /home/mtader/bin/jdftx',
#    executable='/home/mtader/bin/jdftx -c 2',
    #executable='/home/mtader/JDFTx/poly_sbatch.sh',
    #executable='/home/mkelley/jdftx/build/jdftx',
    #executable='/home/mtader/JDFTx/jdftx/',
#    executable='/home/jfm343/JDFTXDIR/build/jdftx',
    #pseudoDir='/home/mtader/JDFTx/jdftx/pseudopotentials/',
    pseudoDir='/home/jfm343/JDFTXDIR/build/pseudopotentials/',
    pseudoSet='GBRV-pbe',
    commands={'elec-cutoff' : '20 100','kpoint': '.5 0 0','kpoint-folding' : '2 2 2', 'electronic-scf':' ', 'elec-smearing': 'Fermi 0.01'}
    #commands={'elec-cutoff' : '20 100'}
)

images = [initial]
for i in range(2):
#for i in range(3):
    image = initial.copy()
    image.set_calculator(JDFTx())
    image.set_constraint(constraint)
    images.append(image)

images.append(final)

neb = NEB(images, parallel=True)
neb.interpolate()
qn = BFGS(neb, logfile='qn.log', trajectory='nebpar_moretraj.traj')
#qn = BFGS(neb, trajectory='neb.traj')
qn.run(fmax=0.05)
