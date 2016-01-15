from proteus.default_so import *
from proteus import Context

name = "tank"

case = __import__(name)
Context.setFromModule(case)
ct = Context.get()

if ct.useOnlyVF:
    pnList = [("twp_navier_stokes_p",   #0
               "twp_navier_stokes_n"),
              ("vof_p",                 #1
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
    
if tank.useRANS > 0:
    pnList += [("kappa_p",
               "kappa_n")]
    pnList += [("dissipation_p",
                "dissipation_n")]
    
# systemStepControllerType = ISO_fixed_MinAdaptiveModelStep
systemStepControllerType = Sequential_MinAdaptiveModelStep

needEBQ_GLOBAL = False
needEBQ = False

tnList = [0.0,ct.dt_init]+[i*ct.dt_fixed for i in range(1,ct.nDTout+1)] 
#archiveFlag = ArchiveFlags.EVERY_SEQUENCE_STEP
