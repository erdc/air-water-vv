from builtins import str
from builtins import range
import sharp_crested_weir
from proteus.default_so import *

if sharp_crested_weir.useOnlyVF:
    pnList = [("twp_navier_stokes_p", "twp_navier_stokes_n"),
              ("vof_p",               "vof_n")]
else:
    pnList = [("twp_navier_stokes_p", "twp_navier_stokes_n"),
              ("vof_p",               "vof_n"),
              ("ls_p",                "ls_n"),
              ("redist_p",            "redist_n"),
              ("ls_consrv_p",         "ls_consrv_n")]
    
    
if sharp_crested_weir.useRANS > 0:
    pnList.append(("kappa_p",
                   "kappa_n"))
    pnList.append(("dissipation_p",
                   "dissipation_n"))
name = "sharp_crested_weir_p" 

if sharp_crested_weir.timeDiscretization == 'flcbdf':
    systemStepControllerType = Sequential_MinFLCBDFModelStep
    systemStepControllerType = Sequential_MinAdaptiveModelStep
else:
    systemStepControllerType = Sequential_MinAdaptiveModelStep

needEBQ_GLOBAL = False
needEBQ = False

tnList = [0.0, sharp_crested_weir.dt_init] + [i * sharp_crested_weir.dt_fixed for i in range(1,
                                                                                             sharp_crested_weir.nDTout + 1)]

info = open("TimeList.txt","w")


for time in tnList:
    info.write(str(time)+"\n")
info.close()
