from ase.constraints import FixAtoms
from ase import Atoms
from JDFTx import JDFTx
from ase.optimize import QuasiNewton

MOH1 = Atoms('HHHHHCO',
    positions=[(0.979411783123183,6.414138935650924,-2.011166291603455),(0.521689756694117,   4.061935011164090,  -0.343251071565552),(-0.558015842556813,   0.595303680455614,   1.881192774202853),(1.936657168198373,  -0.752252269253269,  -0.481032628182144),(-1.029905078220830,  -0.530088128172709,  -1.299792417580344),(0.000000000000000,   0.000000000000000,   0.000000000000000),(0.712178198980380,   5.912323472421656,  -0.251849528512894)],
    cell=[(15,0,0),(0,25,0),(0,0,15)]
)


MOH2 = Atoms('HHHHHCO',
    positions=[(0.080000377422581,   7.990319110792399,  -0.150240606895511),(-0.063678737898889,   2.076648864901587,  -0.132043177143289),(-0.843073352078655, -0.794308467475609,   1.912870110018682),(1.936657168198373,  -0.752252269253269,  -0.481032628182144),(-1.029905078220830,  -0.530088128172709,  -1.299792417580344),(0.000000000000000,   0.000000000000000,   0.000000000000000),(-0.080000377422581,   6.549999056443548,   0.975059915834764)],    cell=[(15,0,0),(0,25,0),(0,0,15)])

c = FixAtoms(mask=[0, 0, 0, 0, 1, 1, 1])
MOH1.set_constraint(c)
MOH2.set_constraint(c)

calculator = JDFTx(
    executable='/home/mtader/bin/jdftx',
    #executable='/home/mkelley/jdftx/build/jdftx',
    #executable='/home/mtader/JDFTx/jdftx/',
#    executable='/home/jfm343/JDFTXDIR/build/jdftx',
    #pseudoDir='/home/mtader/JDFTx/jdftx/pseudopotentials/',
    pseudoDir='/home/jfm343/JDFTXDIR/build/pseudopotentials/',
    pseudoSet='GBRV-pbe',
    commands={'elec-cutoff' : '20 100','kpoint': '.5 0 0','kpoint-folding' : '2 2 2', 'electronic-scf':' ', 'elec-smearing': 'Fermi 0.01'}
    #commands={'elec-cutoff' : '20 100'}
)


MOH1.set_calculator(JDFTx())
MOH2.set_calculator(JDFTx())

qn = QuasiNewton(MOH1, trajectory='initial.traj')
qn.run(fmax=0.05)
qn = QuasiNewton(MOH2, trajectory='final.traj')
qn.run(fmax=0.05)
