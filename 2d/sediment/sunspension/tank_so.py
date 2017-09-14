from proteus.default_so import *
import tank
import os


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



from proteus.SplitOperator import Sequential_FixedStep_Simple, defaultSystem

if tank.sedimentDynamics:
    pnList = [("vos_p",               "vos_n"),#0
              ("vof_p",               "vof_n"),#1
              ("ls_p",                "ls_n"),#2
              ("redist_p",            "redist_n"),#3
              ("ls_consrv_p",         "ls_consrv_n"),#4
              ("threep_navier_stokes_sed_p", "threep_navier_stokes_sed_n"),#5
              ("twp_navier_stokes_p", "twp_navier_stokes_n"),#6
              ("pressureincrement_p", "pressureincrement_n"),#7
              ("pressure_p", "pressure_n"),#8
              ("pressureInitial_p", "pressureInitial_n")]#9
    tank.VOS_model=0
    tank.VOF_model=1
    tank.LS_model=2
    tank.RD_model=3
    tank.MCORR_model=4
    tank.SED_model=5
    tank.V_model=6
    tank.PINC_model=7
    tank.PRESSURE_model=8
    tank.PINIT_model=9
else:
    pnList = [("vof_p",               "vof_n"),#0
              ("ls_p",                "ls_n"),#1
              ("redist_p",            "redist_n"),#2
              ("ls_consrv_p",         "ls_consrv_n"),#3
              ("twp_navier_stokes_p", "twp_navier_stokes_n"),#4
              ("pressureincrement_p", "pressureincrement_n"),#5
              ("pressure_p", "pressure_n"),#6
              ("pressureInitial_p", "pressureInitial_n")]#7
    tank.VOS_model=None
    tank.SED_model=None
    tank.VOF_model=0
    tank.LS_model=1
    tank.RD_model=2
    tank.MCORR_model=3
    tank.V_model=4
    tank.PINC_model=5
    tank.PRESSURE_model=6
    tank.PINIT_model=7

if tank.useRANS > 0:
    pnList.append(("kappa_p",
                   "kappa_n"))
    pnList.append(("dissipation_p",
                   "dissipation_n"))
name = "tank"

#modelSpinUpList = [tank.VOF_model, tank.LS_model, tank.V_model, tank.PINIT_model]
modelSpinUpList = [tank.PINIT_model]

class Sequential_MinAdaptiveModelStepPS(Sequential_MinAdaptiveModelStep):
    def __init__(self,modelList,system=defaultSystem,stepExact=True):
        Sequential_MinAdaptiveModelStep.__init__(self,modelList,system,stepExact)
        self.modelList = modelList[:len(pnList)-1]


systemStepControllerType = Sequential_MinAdaptiveModelStepPS

needEBQ_GLOBAL = False
needEBQ = False

tnList = [0.0,tank.dt_init]+[i*tank.dt_fixed for i in range(1,tank.nDTout+1)]

info = open("TimeList.txt","w")
#archiveFlag = ArchiveFlags.EVERY_SEQUENCE_STEP
