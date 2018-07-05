from builtins import str
from builtins import range
from proteus.default_so import *
import dambreak_Ubbink_fine

if dambreak_Ubbink_fine.useOnlyVF:
    pnList = [("twp_navier_stokes_p", "twp_navier_stokes_n"),
              ("vof_p",               "vof_n")]
else:
    pnList = [("twp_navier_stokes_p", "twp_navier_stokes_n"),
              ("vof_p",               "vof_n"),
              ("ls_p",                "ls_n"),
              ("redist_p",            "redist_n"),
              ("ls_consrv_p",         "ls_consrv_n")]
    
    
if dambreak_Ubbink_fine.useRANS > 0:
    pnList.append(("kappa_p",
                   "kappa_n"))
    pnList.append(("dissipation_p",
                   "dissipation_n"))
name = "dambreak_Ubbink_fine_p" 

if dambreak_Ubbink_fine.timeDiscretization == 'flcbdf':
    systemStepControllerType = Sequential_MinFLCBDFModelStep
    systemStepControllerType = Sequential_MinAdaptiveModelStep
else:
    systemStepControllerType = Sequential_MinAdaptiveModelStep

needEBQ_GLOBAL = False
needEBQ = False

tnList = [0.0,dambreak_Ubbink_fine.dt_init]+[i*dambreak_Ubbink_fine.dt_fixed for i in range(1,dambreak_Ubbink_fine.nDTout+1)] 

info = open("TimeList.txt","w")


for time in tnList:
    info.write(str(time)+"\n")
info.close()
