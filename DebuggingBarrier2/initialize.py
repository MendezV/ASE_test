from ase import Atoms
from JDFTx import JDFTx
from ase.parallel import rank, size
from ase.optimize import BFGS
from ase.constraints import FixAtoms



###################################################
###################################################
####################initial state##################
###################################################
###################################################


d =  1. #Initial bond length guess
CO_ini = Atoms('HH', 
    positions=[(0, 0, 0), (0, 0, d)],
    cell=[(6,0,0),(0,6,0),(0,0,7)],
    pbc = [False, False, False])


maski=[atom.tag for atom in CO_ini]
indi=[atom.index for atom in CO_ini]
print([atom.tag for atom in CO_ini],[atom.symbol for atom in CO_ini],[atom.index for atom in CO_ini])


#Set up JDFTx calculator
calculator = JDFTx(
    executable='/home/jfm343/JDFTXDIR/build/jdftx',
    pseudoDir='/home/jfm343/JDFTXDIR/build/pseudopotentials',
    pseudoSet='GBRV-pbe',
    commands={'elec-cutoff' : '20 100', 'elec-ex-corr' :'gga-PBEsol', 'spintype' : 'z-spin','elec-smearing' :'Fermi 0.01' })
CO_ini.set_calculator(calculator)
CO_ini.set_constraint(FixAtoms(indices=indi))

#Structure optimization
dyn1 = BFGS(CO_ini, trajectory='initial.traj')
dyn1.run(fmax=0.05)

calculator.clean()  #Clean up run files from /tmp







###################################################
###################################################
####################final state####################
###################################################
###################################################

d = 1.2 #Final bond length guess
CO_fin = Atoms('HH', 
    positions=[(0, 0, 0), (0, 0, d)],
    cell=[(6,0,0),(0,6,0),(0,0,7)],
    pbc = [False, False, False])



maski=[atom.tag for atom in CO_fin]
indi=[atom.index for atom in CO_fin]
print([atom.tag for atom in CO_fin],[atom.symbol for atom in CO_fin],[atom.index for atom in CO_fin])



#Set up JDFTx calculator
calculator = JDFTx(
    executable='/home/jfm343/JDFTXDIR/build/jdftx',
    pseudoDir='/home/jfm343/JDFTXDIR/build/pseudopotentials',
    pseudoSet='GBRV-pbe',
    commands={'elec-cutoff' : '20 100', 'elec-ex-corr' :'gga-PBEsol', 'spintype' : 'z-spin','elec-smearing' :'Fermi 0.01' })

CO_fin.set_calculator(calculator)
CO_fin.set_constraint(FixAtoms(indices=indi))

#Structure optimization
dyn2 = BFGS(CO_fin, trajectory='final.traj')
dyn2.run(fmax=0.05)

calculator.clean()  #Clean up run files from /tmp
