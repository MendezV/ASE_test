from ase import Atoms
from JDFTx import JDFTx
from ase.parallel import rank, size
from ase.optimize import BFGS
from ase.constraints import FixAtoms

d = 1.1 #Initial bond length guess
CO = Atoms('CO', 
    positions=[(0, 0, 0), (0, 0, d)],
    cell=[(6,0,0),(0,6,0),(0,0,7)],
    pbc = [False, False, False])



mask = [atom.tag for atom in CO]
#print(mask)
CO.set_constraint(FixAtoms(mask=mask))


#Set up JDFTx calculator
calculator = JDFTx(
    executable='srun -n 1 -N 1 -c 12 --exclude=node[1001-1032] /home/jfm343/JDFTXDIR/build/jdftx',
    pseudoDir='/home/jfm343/JDFTXDIR/build/pseudopotentials',
    pseudoSet='GBRV-pbe',
    commands={'elec-cutoff' : '20 100', 'elec-ex-corr' :'gga-PBEsol'})
CO.set_calculator(calculator)

#Structure optimization
dyn = BFGS(CO, trajectory='initial.traj')
dyn.run(fmax=0.05)

calculator.clean()  #Clean up run files from /tmp

