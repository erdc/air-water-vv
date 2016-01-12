from proteus.default_so import *
from proteus import Context

name = "tank"

case = __import__(name)
Context.setFromModule(case)
ct = Context.get()

if ct.useOnlyVF:
    pnList = [("twp_navier_stokes_p",  #0
               "twp_navier_stokes_n"),
              ("vof_p",                #1
               "vof_n")]
else:
    pnList = [("twp_navier_stokes_p",  #0
               "twp_navier_stokes_n"),
              ("vof_p",                #1
               "vof_n"),
              ("ls_p",                 #2
               "ls_n"),
              ("redist_p",             #3
               "redist_n"),
              ("ls_consrv_p",          #4
               "ls_consrv_n")]
    
if ct.movingDomain:
    pnList = [("moveMesh_p","moveMesh_n")]+pnList
    
if ct.useRANS > 0:
    pnList += [("kappa_p",
                "kappa_n")]
    pnList += [("dissipation_p",
                "dissipation_n")] 


if ct.timeDiscretization == 'flcbdf':
    systemStepControllerType = Sequential_MinFLCBDFModelStep
    systemStepControllerType = Sequential_MinAdaptiveModelStep
else:
    systemStepControllerType = Sequential_MinAdaptiveModelStep      # systemStepControllerType = ISO_fixed_MinAdaptiveModelStep


needEBQ_GLOBAL = False
needEBQ = False


tnList = [0.0,0.001,1] #ct.dt_init]+[i*ct.dt_fixed for i in range(1,ct.nDTout+1)] 

