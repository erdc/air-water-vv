from proteus.default_so import *
import wavesloshing_laminar_unstruct_coarse

if wavesloshing_laminar_unstruct_coarse.useOnlyVF:
    pnList = [("twp_navier_stokes_p", "twp_navier_stokes_n"),
              ("vof_p",               "vof_n")]
else:
    pnList = [("twp_navier_stokes_p", "twp_navier_stokes_n"),
              ("vof_p",               "vof_n"),
              ("ls_p",                "ls_n"),
              ("redist_p",            "redist_n"),
              ("ls_consrv_p",         "ls_consrv_n")]
    
    
if wavesloshing_laminar_unstruct_coarse.useRANS > 0:
    pnList.append(("kappa_p",
                   "kappa_n"))
    pnList.append(("dissipation_p",
                   "dissipation_n"))
name = "wavesloshing_laminar_unstruct_coarse_p" 

if wavesloshing_laminar_unstruct_coarse.timeDiscretization == 'flcbdf':
    systemStepControllerType = Sequential_MinFLCBDFModelStep
    systemStepControllerType = Sequential_MinAdaptiveModelStep
else:
    systemStepControllerType = Sequential_MinAdaptiveModelStep

needEBQ_GLOBAL = False
needEBQ = False

tnList = [0.0,wavesloshing_laminar_unstruct_coarse.dt_init]+[i*wavesloshing_laminar_unstruct_coarse.dt_fixed for i in range(1,wavesloshing_laminar_unstruct_coarse.nDTout+1)] 

info = open("TimeList.txt","w")


