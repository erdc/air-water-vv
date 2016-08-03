from proteus.default_so import *
import Mase_and_Kirby_shoaling

if Mase_and_Kirby_shoaling.useOnlyVF:
    pnList = [("twp_navier_stokes_p", "twp_navier_stokes_n"),
              ("vof_p",               "vof_n")]
else:
    pnList = [("twp_navier_stokes_p", "twp_navier_stokes_n"),
              ("vof_p",               "vof_n"),
              ("ls_p",                "ls_n"),
              ("redist_p",            "redist_n"),
              ("ls_consrv_p",         "ls_consrv_n")]
    
    
if Mase_and_Kirby_shoaling.useRANS > 0:
    pnList.append(("kappa_p",
                   "kappa_n"))
    pnList.append(("dissipation_p",
                   "dissipation_n"))
name = "tank_p" 

if Mase_and_Kirby_shoaling.timeDiscretization == 'flcbdf':
    systemStepControllerType = Sequential_MinFLCBDFModelStep
    systemStepControllerType = Sequential_MinAdaptiveModelStep
else:
    systemStepControllerType = Sequential_MinAdaptiveModelStep

needEBQ_GLOBAL = False
needEBQ = False

tnList = [0.0, Mase_and_Kirby_shoaling.dt_init] + [i * Mase_and_Kirby_shoaling.dt_fixed for i in range(1, Mase_and_Kirby_shoaling.nDTout + 1)]
archiveFlag = ArchiveFlags.EVERY_SEQUENCE_STEP
