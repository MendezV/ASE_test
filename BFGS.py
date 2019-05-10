from ase import Atoms
from JDFTx import JDFTx
from ase.parallel import rank, size

d = 1.1 #Initial bond length guess
CO = Atoms('CO', 
    positions=[(0, 0, 0), (0, 0, d)],
    cell=[(6,0,0),(0,6,0),(0,0,7)],
    pbc = [False, False, False])

#Set up JDFTx calculator
calculator = JDFTx(
    executable='srun -n 1 -N 1 -c 12 --exclude=node[1001-1032] /home/jfm343/JDFTXDIR/build/jdftx',
    pseudoDir='/home/jfm343/JDFTXDIR/build/pseudopotentials',
    pseudoSet='GBRV-pbe',
    commands={'elec-cutoff' : '20 100', 'elec-ex-corr' :'gga-PBEsol'})
CO.set_calculator(calculator)

#Structure optimization
from ase.optimize import BFGS
dyn = BFGS(CO)
dyn.run(fmax=0.05)

calculator.clean()  #Clean up run files from /tmp
