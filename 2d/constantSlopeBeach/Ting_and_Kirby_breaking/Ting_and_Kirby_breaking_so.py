from proteus.default_so import *
import Ting_and_Kirby_breaking

if Ting_and_Kirby_breaking.useOnlyVF:
    pnList = [("twp_navier_stokes_p", "twp_navier_stokes_n"),
              ("vof_p",               "vof_n")]
else:
    pnList = [("twp_navier_stokes_p", "twp_navier_stokes_n"),
              ("vof_p",               "vof_n"),
              ("ls_p",                "ls_n"),
              ("redist_p",            "redist_n"),
              ("ls_consrv_p",         "ls_consrv_n")]
    
    
if Ting_and_Kirby_breaking.useRANS > 0:
    pnList.append(("kappa_p",
                   "kappa_n"))
    pnList.append(("dissipation_p",
                   "dissipation_n"))
name = "tank_p" 

if Ting_and_Kirby_breaking.timeDiscretization == 'flcbdf':
    systemStepControllerType = Sequential_MinFLCBDFModelStep
    systemStepControllerType = Sequential_MinAdaptiveModelStep
else:
    systemStepControllerType = Sequential_MinAdaptiveModelStep

needEBQ_GLOBAL = False
needEBQ = False

tnList = [0.0, Ting_and_Kirby_breaking.dt_init] + [i * Ting_and_Kirby_breaking.dt_fixed for i in range(1, Ting_and_Kirby_breaking.nDTout + 1)]
archiveFlag = ArchiveFlags.EVERY_SEQUENCE_STEP
