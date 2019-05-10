from JDFTx import JDFTx
from ase import Atoms
from ase.io import read
from ase.constraints import FixAtoms
from ase.neb import NEB
from ase.optimize import BFGS,LBFGS

initial = Atoms('HHHHHCO',positions=[(0.979411783123183,6.414138935650924,-2.011166291603455),(0.521689756694117,   4.061935011164090,  -0.343251071565552),(-0.558015842556813,   0.595303680455614,   1.881192774202853),(1.936657168198373,  -0.752252269253269,  -0.481032628182144),(-1.029905078220830,  -0.530088128172709,  -1.299792417580344),(0.000000000000000,   0.000000000000000,   0.000000000000000),(0.712178198980380,   5.912323472421656,  -0.251849528512894)],cell=[(15,0,0),(0,25,0),(0,0,15)])
final = Atoms('HHHHHCO',positions=[(0.080000377422581,   7.990319110792399,  -0.150240606895511),(-0.063678737898889,   2.076648864901587,  -0.132043177143289),(-0.843073352078655, -0.794308467475609,   1.912870110018682),(1.936657168198373,  -0.752252269253269,  -0.481032628182144),(-1.029905078220830,  -0.530088128172709,  -1.299792417580344),(0.000000000000000,   0.000000000000000,   0.000000000000000),(-0.080000377422581,   6.549999056443548,   0.975059915834764)],    cell=[(15,0,0),(0,25,0),(0,0,15)])

c2 = FixAtoms(indices=[4,5,6])
images = [initial]
for i in range(3):
    image = initial.copy()
    image.set_calculator(JDFTx(executable='/home/jfm343/JDFTXDIR/build/jdftx',pseudoDir='/home/jfm343/JDFTXDIR/build/pseudopotentials',pseudoSet='GBRV-pbe',commands={'elec-cutoff' : '20 100','kpoint': '.5 0 0','kpoint-folding' : '2 2 2', 'electronic-scf':' ', 'elec-smearing': 'Fermi 0.01'}))
    images.append(image)
    image.set_constraint(c2) ## so that beads are able to move in intermediate steps

images.append(final)

neb = NEB(images, parallel=True)
neb.interpolate()
qn = LBFGS(neb, trajectory='neb.traj')
qn.run(fmax=0.05)

