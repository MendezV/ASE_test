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


ini = Atoms('HHHHHCO',positions=[(0.080000377422581,   7.990319110792399,  -0.150240606895511),(-0.063678737898889,   2.076648864901587,  -0.132043177143289),(-0.843073352078655, -0.794308467475609,   1.912870110018682),(1.936657168198373,  -0.752252269253269,  -0.481032628182144),(-1.029905078220830,  -0.530088128172709,  -1.299792417580344),(0.000000000000000,   0.000000000000000,   0.000000000000000),(-0.080000377422581,   6.549999056443548,   0.975059915834764)],    cell=[(7.4084,0,0),(0,9.525,0),(0,0,7.4084)],pbc = [True, True, True])

maski=[atom.tag for atom in ini]
indi=[atom.index for atom in ini]
print([atom.tag for atom in ini],[atom.symbol for atom in ini],[atom.index for atom in ini])
c = FixAtoms(mask=[0, 0, 0, 0, 1, 1, 1])
c2 = FixAtoms(indices=[3,4,5])

#Set up JDFTx calculator
calculator = JDFTx(
    executable='srun -n 1 -N 1 -c 12 --exclude=node[1001-1032] /home/jfm343/JDFTXDIR/build/jdftx',
    pseudoDir='/home/jfm343/JDFTXDIR/build/pseudopotentials',
    pseudoSet='GBRV-pbe',
    commands={'elec-cutoff' : '20 100','kpoint': '.5 0 0','kpoint-folding' : '4 1 1', 'electronic-minimize':'nIterations 100 ', 'elec-smearing': 'Fermi 0.001','elec-n-bands': '23','spintype' :'z-spin','ionic-minimize' : 'nIterations 10 ' })
ini.set_calculator(calculator)
ini.set_constraint(c2)

#Structure optimization
dyn1 = BFGS(ini, trajectory='initial.traj')
dyn1.run(fmax=0.05)

calculator.clean()  #Clean up run files from /tmp







###################################################
###################################################
####################final state####################
###################################################
###################################################


fin = Atoms('HHHHHCO',positions=[(0.979411783123183,6.414138935650924,-2.011166291603455),(0.521689756694117,   4.061935011164090,  -0.343251071565552),(-0.558015842556813,   0.595303680455614,   1.881192774202853),(1.936657168198373,  -0.752252269253269,  -0.481032628182144),(-1.029905078220830,  -0.530088128172709,  -1.299792417580344),(0.000000000000000,   0.000000000000000,   0.000000000000000),(0.712178198980380,   5.912323472421656,  -0.251849528512894)],cell=[(7.4084,0,0),(0,9.525,0),(0,0,7.4084)],pbc = [True, True, True])


maski=[atom.tag for atom in fin]
indi=[atom.index for atom in fin]
print([atom.tag for atom in fin],[atom.symbol for atom in fin],[atom.index for atom in fin])
c = FixAtoms(mask=[0, 0, 0, 0, 1, 1, 1])
c2 = FixAtoms(indices=[3,4,5])

#Set up JDFTx calculator
calculator = JDFTx(
    executable='srun -n 1 -N 1 -c 12 --exclude=node[1001-1032] /home/jfm343/JDFTXDIR/build/jdftx',
    pseudoDir='/home/jfm343/JDFTXDIR/build/pseudopotentials',
    pseudoSet='GBRV-pbe',
    commands={'elec-cutoff' : '20 100','kpoint': '.5 0 0','kpoint-folding' : '4 1 1', 'electronic-scf':' ', 'elec-smearing': 'Fermi 0.001','elec-n-bands': '23','spintype' :'z-spin'})

fin.set_calculator(calculator)
fin.set_constraint(c2)

#Structure optimization
dyn2 = BFGS(fin, trajectory='final.traj')
dyn2.run(fmax=0.05)

calculator.clean()  #Clean up run files from /tmp

