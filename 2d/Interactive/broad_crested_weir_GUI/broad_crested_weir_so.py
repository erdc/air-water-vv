"""
Split operator module for two-phase flow
"""
import os
from proteus.default_so import *
from proteus import Context


# Create context from main module
name_so = os.path.basename(__file__)
if '_so.py' in name_so[-6:]:
    name = name_so[:-6]
elif '_so.pyc' in name_so[-7:]:
    name = name_so[:-7]
else:
    raise NameError, 'Split operator module must end with "_so.py"'

try:
    case = __import__(name)
    Context.setFromModule(case)
    ct = Context.get()
except ImportError:
    raise ImportError, str(name) + '.py not found'

if ct.useOnlyVF:
    pnList = [("twp_navier_stokes_p", "twp_navier_stokes_n"),
              ("vof_p",               "vof_n")]
else:
    pnList = [("twp_navier_stokes_p", "twp_navier_stokes_n"),
              ("vof_p",               "vof_n"),
              ("ls_p",                "ls_n"),
              ("redist_p",            "redist_n"),
              ("ls_consrv_p",         "ls_consrv_n")]

if ct.useRANS > 0:
    pnList.append(("kappa_p",
                   "kappa_n"))
    pnList.append(("dissipation_p",
                   "dissipation_n"))
name = "broad_crested_weir_p" 

if ct.dt_fixed:
#    systemStepControllerType = Sequential_FixedStep
    systemStepControllerType = Sequential_MinAdaptiveModelStep
    dt_system_fixed = ct.dt_fixed
    stepExactSystem=False
else:  # use CFL
    systemStepControllerType = Sequential_MinAdaptiveModelStep
    stepExactSystem=False

needEBQ_GLOBAL = False
needEBQ = False

tnList = [0.0]+[i*ct.dt_fixed for i in range(1,ct.nDTout+1)]

#Disabled to allow GUI to create images without restarting simulation
info = open("TimeList.txt","w")

for time in tnList:
   info.write(str(time)+"\n")
info.close()
