from proteus import *
from proteus.default_p import *
from broad_crested_weir import *
from proteus.mprans import RANS2P
import BC as BC 
setBC = BC.boundaryConditions()
def createBoundaryCondition(x,flag,BCType):
    if flag == boundaryTags['top']:
        return setBC.hydrostaticPressureOutletWithDepth(BCType,x,waterLine_z,rho_1,rho_0,g,refLevel=L[1],b_or=[0,1,0])
    if flag == boundaryTags['left']:
        return setBC.twoPhaseVelocityInlet(BCType,x,
                                     U=[inflow_velocity,0,0],
                                     seaLevel=waterLine_z,
                                     b_or=[1,0,0]
                                     )
                                     
    if flag == boundaryTags['right']:
        return setBC.hydrostaticPressureOutletWithDepth(BCType,x,waterLine_z,rho_1,rho_0,g,refLevel=L[1],b_or=[1,0,0])
    if flag == boundaryTags['bottom']:
        return setBC.freeSlip(BCType)
