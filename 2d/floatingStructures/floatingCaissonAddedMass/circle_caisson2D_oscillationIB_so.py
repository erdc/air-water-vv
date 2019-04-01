"""
Split operator module for two-phase flow
"""

import os
from proteus.default_so import *
from proteus import Context

# Create context from main module
name_so = os.path.basename(__file__)
if '_so.py' in name_so[-6:]:
    name = name_so[:-6]
elif '_so.pyc' in name_so[-7:]:
    name = name_so[:-7]
else:
    raise(NameError, 'Split operator module must end with "_so.py"')

case = __import__(name)
Context.setFromModule(case)
ct = Context.get()

# List of p/n files
pnList = [("vof3p_p", "vof3p_n"),#0
          ("ls3p_p", "ls3p_n"),#1
          ("redist3p_p", "redist3p_n"),#2
          ("ls_consrv3p_p", "ls_consrv3p_n"),#3
          ("rans3p_p", "rans3p_n"),#4
          ("pressureincrement_p", "pressureincrement_n"),#5
          ("pressure_p", "pressure_n"),#6
          ("added_massIB_p","added_massIB_n"),#7
          ("pressureInitial_p", "pressureInitial_n")]#8
modelSpinUpList = [8]
#systemStepControllerType = ISO_fixed_MinAdaptiveModelStep
#if ct.dt_fixed:
#    systemStepControllerType = Sequential_FixedStep
#    dt_system_fixed = ct.dt_fixed
#    systemStepExact=False
#else:  # use CFL
#    systemStepControllerType = Sequential_MinAdaptiveModelStep
#    systemStepExact=False
#class Sequential_PS(Sequential_FixedStep):
#    def __init__(self,modelList,system=defaultSystem,stepExact=True):
#        Sequential_FixedStep.__init__(self,modelList,system,stepExact)
#        self.modelList = modelList[:len(pnList)]

class Sequential_FixedStepPS(Sequential_FixedStep):
    def __init__(self,modelList,system=defaultSystem,stepExact=True):
        Sequential_FixedStep.__init__(self,modelList,system,stepExact)
        self.modelList = modelList[:len(pnList)-1]
systemStepControllerType = Sequential_FixedStepPS
systemStepExact=False
dt_system_fixed = ct.opts.dt_fixed

needEBQ_GLOBAL = False
needEBQ = False

#modelSpinUpList = [0]  # for initial conditions of movemesh
#archiveFlag = ArchiveFlags.EVERY_SEQUENCE_STEP

if ct.opts.nsave == 0:
    if ct.dt_fixed > 0:
        archiveFlag = ArchiveFlags.EVERY_USER_STEP
        if ct.dt_init < ct.dt_fixed:
            tnList = [0., ct.dt_init, ct.dt_fixed, ct.T]
        else:
            tnList = [0., ct.dt_fixed, ct.T]
    else:
          tnList = [0., ct.dt_init, ct.T]
else:
    tnList=[0.0,ct.dt_init]+[ct.dt_init+ i*ct.dt_out for i in range(1,ct.nDTout+1)]

