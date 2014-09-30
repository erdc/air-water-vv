from proteus.default_so import *
import dambreak_Colagrossi_medium_R_40

if dambreak_Colagrossi_medium_R_40.useOnlyVF:
    pnList = [("twp_navier_stokes_p", "twp_navier_stokes_n"),
              ("vof_p",               "vof_n")]
else:
    pnList = [("twp_navier_stokes_p", "twp_navier_stokes_n"),
              ("vof_p",               "vof_n"),
              ("ls_p",                "ls_n"),
              ("redist_p",            "redist_n"),
              ("ls_consrv_p",         "ls_consrv_n")]
    
    
if dambreak_Colagrossi_medium_R_40.useRANS > 0:
    pnList.append(("kappa_p",
                   "kappa_n"))
    pnList.append(("dissipation_p",
                   "dissipation_n"))
name = "dambreak_Colagrossi_medium_R_40_p" 

if dambreak_Colagrossi_medium_R_40.timeDiscretization == 'flcbdf':
    systemStepControllerType = Sequential_MinFLCBDFModelStep
    systemStepControllerType = Sequential_MinAdaptiveModelStep
else:
    systemStepControllerType = Sequential_MinAdaptiveModelStep

needEBQ_GLOBAL = False
needEBQ = False

tnList = [0.0,dambreak_Colagrossi_medium_R_40.dt_init]+[i*dambreak_Colagrossi_medium_R_40.dt_fixed for i in range(1,dambreak_Colagrossi_medium_R_40.nDTout+1)] 

info = open("TimeList.txt","w")


