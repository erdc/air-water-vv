from proteus.default_so import *
from proteus import Context
import tank
Context.setFromModule(tank)
ct = Context.get()

if ct.useOnlyVF:
    pnList = [("twp_navier_stokes_p", "twp_navier_stokes_n"),
              ("vof_p",               "vof_n")]
else:
    pnList = [("twp_navier_stokes_p", "twp_navier_stokes_n"),
              ("vof_p",               "vof_n"),
              ("ls_p",                "ls_n"),
              ("redist_p",            "redist_n"),
              ("ls_consrv_p",         "ls_consrv_n")]
    
    
if ct.useRANS > 0:
    pnList.append(("kappa_p",
                   "kappa_n"))
    pnList.append(("dissipation_p",
                   "dissipation_n"))
name = "tank_p" 

if ct.timeDiscretization == 'flcbdf':
    systemStepControllerType = Sequential_MinFLCBDFModelStep
    systemStepControllerType = Sequential_MinAdaptiveModelStep
else:
    systemStepControllerType = Sequential_MinAdaptiveModelStep

needEBQ_GLOBAL = False
needEBQ = False

tnList = [0.0,ct.dt_init]+[i*ct.dt_fixed for i in range(1,ct.nDTout+1)] 
