from proteus.default_so import *
import dambreak


if dambreak.ct.only_vof:
    pnList = [
                ("twp_navier_stokes_p", "twp_navier_stokes_n"),
                ("vof_p",               "vof_n"),
            ]
    dambreak.V_model=0
    dambreak.VF_model=1
    dambreak.LS_model=None
    dambreak.RD_model=None
    dambreak.MCORR_model=None
else:
    pnList = [
                ("twp_navier_stokes_p", "twp_navier_stokes_n"),
                ("vof_p",               "vof_n"),
                ("ls_p",                "ls_n"),
                ("redist_p",            "redist_n"),
                ("ls_consrv_p",         "ls_consrv_n"),
            ]
    dambreak.V_model=0
    dambreak.VF_model=1
    dambreak.LS_model=2
    dambreak.RD_model=3
    dambreak.MC_model=4


# if dambreak.useRANS > 0:
#     pnList.append(("kappa_p",
#                    "kappa_n"))
#     pnList.append(("dissipation_p",
                #    "dissipation_n"))
name = "dambreak_p" 

if dambreak.timeDiscretization == 'flcbdf':
    systemStepControllerType = Sequential_MinFLCBDFModelStep
    systemStepControllerType = Sequential_MinAdaptiveModelStep
else:
    systemStepControllerType = Sequential_MinAdaptiveModelStep

needEBQ_GLOBAL = False
needEBQ = False

tnList = [0.0,dambreak.dt_init]+[i*dambreak.dt_fixed for i in range(1,dambreak.nDTout+1)] 
archiveFlag = ArchiveFlags.EVERY_USER_STEP
