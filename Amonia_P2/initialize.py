from ase import Atoms
from JDFTx import JDFTx
from ase.parallel import rank, size
from ase.optimize import BFGS
from ase.constraints import FixAtoms
from ase.io import read, write



###################################################
###################################################
####################initial state##################
###################################################
###################################################


Amo_ini = Atoms('NHHH',positions=[(0.257,  -0.363 ,  0.000), (0.257,   0.727 ,  0.000),(0.771 , -0.727  , 0.890),(0.771,  -0.727 , -0.890)],cell=[(8,0,0),(0,8,0),(0,0,8)], pbc = [False, False, False])

maski=[atom.tag for atom in Amo_ini]
indi=[atom.index for atom in Amo_ini]
print([atom.tag for atom in Amo_ini],[atom.symbol for atom in Amo_ini],[atom.index for atom in Amo_ini])


#Set up JDFTx calculator
calculator = JDFTx( executable='srun -n 1 -N 1 -c 12 --exclude=node[1001-1032]  /home/jfm343/JDFTXDIR/build/jdftx',pseudoDir='/home/jfm343/JDFTXDIR/build/pseudopotentials',pseudoSet='GBRV-pbe',commands={'elec-cutoff' : '20 100', 'elec-ex-corr' :'gga-PBEsol', 'spintype' : 'z-spin','elec-smearing' :'Fermi 0.01' })

Amo_ini.set_calculator(calculator)
Amo_ini.set_constraint(FixAtoms(indices=[0]))

#Structure optimization
dyn1 = BFGS(Amo_ini, trajectory='initial.traj')
dyn1.run(fmax=0.05)

calculator.clean()  #Clean up run files from /tmp



###################################################
###################################################
####################final state####################
###################################################
###################################################

d = 1.2 #Final bond length guess
Amo_fini = Atoms('NHHH',positions=[(0.257,  -0.363 ,  0.000), (0.257,   0.727 ,  0.000),(-0.771 , -0.727  , 0.890),(-0.771,  -0.727 , -0.890)],cell=[(8,0,0),(0,8,0),(0,0,8)],pbc = [False, False, False])



maski=[atom.tag for atom in Amo_fini]
indi=[atom.index for atom in Amo_fini]
print([atom.tag for atom in Amo_fini],[atom.symbol for atom in Amo_fini],[atom.index for atom in Amo_fini])



#Set up JDFTx calculator
calculator = JDFTx( executable='srun -n 1 -N 1 -c 12 --exclude=node[1001-1032]  /home/jfm343/JDFTXDIR/build/jdftx',pseudoDir='/home/jfm343/JDFTXDIR/build/pseudopotentials',pseudoSet='GBRV-pbe',commands={'elec-cutoff' : '20 100', 'elec-ex-corr' :'gga-PBEsol', 'spintype' : 'z-spin','elec-smearing' :'Fermi 0.01' })

Amo_fini.set_calculator(calculator)
Amo_fini.set_constraint(FixAtoms(indices=[0]))

#Structure optimization
dyn2 = BFGS(Amo_fini, trajectory='final.traj')
dyn2.run(fmax=0.05)

calculator.clean()  #Clean up run files from /tmp
