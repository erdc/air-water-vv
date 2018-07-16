from proteus.default_p import *
from proteus.mprans import MoveMeshMonitor
import numpy as np
from proteus import Context

ct = Context.get()
domain = ct.domain
nd = ct.domain.nd
mesh = domain.MeshOptions


genMesh = mesh.genMesh
movingDomain = ct.movingDomain
T = ct.T

initialConditions = None

analyticalSolution = {}

#LevelModelType = MoveMesh.LevelModel
coefficients = MoveMeshMonitor.Coefficients(ct.my_func,
                                            nd=ct.domain.nd,
                                            he_max=ct.he_max,
                                            he_min=ct.he_min,
                                            LS_MODEL=3,
                                            ME_MODEL=0,
                                            fixedNodeMaterialTypes=ct.fixedNodes,
                                            nSmoothIn=ct.nSmoothIn,
                                            nSmoothOut=ct.nSmoothOut,
                                            epsTimeStep=ct.epsTimeStep,
                                            epsFact_density=ct.epsFact_density,
                                            grading=ct.grading_mesh,
                                            grading_type=ct.grading_type
                                            )

def getDBC(x,flag):
    return None

dirichletConditions = {0:getDBC}

def getDFBC(x, flag):
    if flag != 0:
        return lambda x, t: 0.

diffusiveFluxBoundaryConditions = {0:{0:getDFBC}}
