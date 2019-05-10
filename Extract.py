from ase import Atoms
from JDFTx import JDFTx, shell

d = 1.1 #Initial bond length guess
CO = Atoms('CO', 
    positions=[(0, 0, 0), (0, 0, d)],
    cell=[(6,0,0),(0,6,0),(0,0,7)],
    pbc = [False, False, False])

#Set up JDFTx calculator
calculator = JDFTx(
    pseudoSet='SG15',
    commands={
        'elec-cutoff' : '30',
        'dump' : 'End ElecDensity'
    }
)
CO.set_calculator(calculator)

#Structure optimization
from ase.optimize import BFGS
dyn = BFGS(CO)
dyn.run(fmax=0.05)

#Create XSF and visualize:
shell('createXSF %s/out CO.xsf %s/n' % ((calculator.runDir,)*2))
shell('VESTA CO.xsf &')
calculator.clean()         #Clean up run files from /tmp
