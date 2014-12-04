from proteus.default_so import *
import quiescent_water_test_gauges

if quiescent_water_test_gauges.useOnlyVF:
    pnList = [("twp_navier_stokes_p", "twp_navier_stokes_n"),
              ("vof_p",               "vof_n")]
else:
    pnList = [("twp_navier_stokes_p", "twp_navier_stokes_n"),
              ("vof_p",               "vof_n"),
              ("ls_p",                "ls_n"),
              ("redist_p",            "redist_n"),
              ("ls_consrv_p",         "ls_consrv_n")]
    
    
if quiescent_water_test_gauges.useRANS > 0:
    pnList.append(("kappa_p",
                   "kappa_n"))
    pnList.append(("dissipation_p",
                   "dissipation_n"))
name = "quiescent_water_test_gauges_p" 

if quiescent_water_test_gauges.timeDiscretization == 'flcbdf':
    systemStepControllerType = Sequential_MinFLCBDFModelStep
    systemStepControllerType = Sequential_MinAdaptiveModelStep
else:
    systemStepControllerType = Sequential_MinAdaptiveModelStep

needEBQ_GLOBAL = False
needEBQ = False

tnList = [0.0,quiescent_water_test_gauges.dt_init]+[i*quiescent_water_test_gauges.dt_fixed for i in range(1,quiescent_water_test_gauges.nDTout+1)] 

info = open("TimeList.txt","w")


