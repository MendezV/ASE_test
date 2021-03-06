from JDFTx import JDFTx
from ase.io import read
from ase.constraints import FixAtoms
from ase.neb import NEB
from ase.optimize import BFGS,LBFGS
import sys

initial = read('initial.traj')
final = read('final.traj')


images = [initial]
for i in range(5):
    image = initial.copy()
    image.set_calculator(JDFTx( executable='srun -n 1 -N 1 -c 12 --exclude=node[1001-1032]  /home/jfm343/JDFTXDIR/build/jdftx',pseudoDir='/home/jfm343/JDFTXDIR/build/pseudopotentials',pseudoSet='GBRV-pbe',commands={'elec-cutoff' : '20 100', 'elec-ex-corr' :'gga-PBEsol', 'spintype' : 'z-spin','elec-smearing' :'Fermi 0.01' ,'core-overlap-check' :'None'}))
    images.append(image)
    image.set_constraint(FixAtoms(indices=[0])) ## so that beads are able to move in intermediate steps

images.append(final)

fmaxi=float(sys.argv[1])
kspr=float(sys.argv[2])

neb = NEB(images,k=kspr,parallel=True)
neb.interpolate('idpp')
qn = BFGS(neb, trajectory='neb_f'+sys.argv[1]+'_k'+sys.argv[2]+'.traj')
qn.run(fmax=fmaxi)
