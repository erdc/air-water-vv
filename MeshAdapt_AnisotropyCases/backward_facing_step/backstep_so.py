from proteus.default_so import *
import backstep
from proteus.MeshAdaptPUMI import MeshAdaptPUMI

if backstep.useOnlyVF:
    pnList = [("twp_navier_stokes_p", "twp_navier_stokes_n"),
              ("vof_p",               "vof_n")]
else:
    pnList = [("twp_navier_stokes_p", "twp_navier_stokes_n"),
              ("vof_p",               "vof_n"),
              ("ls_p",                "ls_n"),
              ("redist_p",            "redist_n"),
              ("ls_consrv_p",         "ls_consrv_n")]


if backstep.useRANS > 0:
    pnList.append(("kappa_p",
                   "kappa_n"))
    pnList.append(("dissipation_p",
                   "dissipation_n"))
name = "backstep_p"

systemStepControllerType = Sequential_FixedStep_Simple
#systemStepControllerType = Sequential_MinAdaptiveModelStep

needEBQ_GLOBAL = False
needEBQ = False

tnList = [0.0,backstep.dt_init,2.0*backstep.dt_init]+[i*backstep.dt_fixed for i in range(1,backstep.nDTout+1)]
