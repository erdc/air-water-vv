from __future__ import absolute_import
from builtins import object
from proteus.default_p import *
from proteus import Context
from proteus.mprans import CLSVOF

# ****************************************** #
# ********** READ FROM setup file ********** #
# ****************************************** #
ct = Context.get()
clsvof_parameters   = ct.clsvof_parameters
initialConditions   = ct.initialConditions
boundaryConditions  = ct.boundaryConditions
nd = ct.nd
clsvof_parameters = ct.clsvof_parameters

# DOMAIN #
domain = ct.domain

# ******************************** #
# ********** PARAMETERS ********** #
# ******************************** #
useMetrics = clsvof_parameters['useMetrics']
epsFactHeaviside = clsvof_parameters['epsFactHeaviside']
epsFactDiract = clsvof_parameters['epsFactDirac']
epsFactRedist = clsvof_parameters['epsFactRedist']
lambdaFact = clsvof_parameters['lambdaFact']
outputQuantDOFs = clsvof_parameters['outputQuantDOFs']
computeMetrics = clsvof_parameters['computeMetrics']
disc_ICs = clsvof_parameters['disc_ICs']

# ************************************ #
# ********** MODEL INDEXING ********** #
# ************************************ #
VOS_model=0
CLSVOF_model=1
V_model=3

# ********************************** #
# ********** COEFFICIENTS ********** #
# ********************************** #
LevelModelType = CLSVOF.LevelModel
coefficients = CLSVOF.Coefficients(VOS_model=VOS_model,
                                   V_model=V_model,
                                   ME_model=CLSVOF_model,
                                   useMetrics=useMetrics,
                                   epsFactHeaviside=epsFactHeaviside,
                                   epsFactDirac=epsFactHeaviside,
                                   epsFactRedist=epsFactRedist,
                                   lambdaFact=lambdaFact,
                                   outputQuantDOFs=outputQuantDOFs,
                                   computeMetrics=computeMetrics,
                                   disc_ICs=disc_ICs)
coefficients.variableNames=['phi']
name="clsvof"

# **************************************** #
# ********** INITIAL CONDITIONS ********** #
# **************************************** #
initialConditions  = {0: initialConditions['clsvof']}

# ***************************************** #    
# ********** BOUNDARY CONDITIONS ********** #
# ***************************************** #
dirichletConditions = {0: boundaryConditions['clsvof_DBC']}
advectiveFluxBoundaryConditions = {0: boundaryConditions['clsvof_AFBC']}
diffusiveFluxBoundaryConditions = {0:{0: boundaryConditions['clsvof_DFBC']}}
