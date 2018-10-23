from proteus.default_so import *
import nxnball
from proteus.SplitOperator import Sequential_FixedStep_Simple, defaultSystem

pnList = [("twp_navier_stokes_p", "twp_navier_stokes_n"),]

name = "onefallingball"

class Sequential_MinAdaptiveModelStepPS(Sequential_FixedStep):
    def __init__(self,modelList,system=defaultSystem,stepExact=True):
        Sequential_FixedStep.__init__(self,modelList,system,stepExact)
        self.modelList = modelList[:len(pnList)]

systemStepControllerType = Sequential_MinAdaptiveModelStepPS

needEBQ_GLOBAL = False
needEBQ = False

tnList = nxnball.tnList
dt_system_fixed = nxnball.dt_fixed

systemStepExact = True
#archiveFlag = ArchiveFlags.EVERY_SEQUENCE_STEP


