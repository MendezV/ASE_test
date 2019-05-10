from JDFTx import JDFTx
from ase.io import read
from ase.constraints import FixAtoms
from ase.neb import NEB
from ase.optimize import BFGS,LBFGS

initial = read('initial.traj')
final = read('final.traj')


images = [initial]
for i in range(5):
    image = initial.copy()
    image.set_calculator(JDFTx(executable='jdftx',pseudoDir='/Users/juanmendezvalderrama/Documents/JDFTX/jdftx-1.4.2/jdftx/pseudopotentials',pseudoSet='GBRV-pbe',commands={'elec-cutoff' : '20 100', 'elec-ex-corr' :'gga-PBEsol', 'spintype' : 'z-spin','elec-smearing' :'Fermi 0.01','core-overlap-check' :'None' }))
    images.append(image)
    image.set_constraint(FixAtoms(indices=[0])) ## so that beads are able to move in intermediate steps

images.append(final)

neb = NEB(images,k=0.05,parallel=True)
neb.interpolate('idpp')
qn = BFGS(neb, trajectory='neb4.traj')
qn.run(fmax=0.1)
