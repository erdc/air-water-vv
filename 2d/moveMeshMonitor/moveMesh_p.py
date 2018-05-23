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
                                            boundaryNormals=ct.boundaryNormals,
                                            nSmooth=ct.nSmooth)


def getDBC(x,flag):
    if flag != 0:
        return None

dirichletConditions = {0:getDBC}

def getDFBC(x, flag):
    if flag != 0:
        return lambda x, t: 0.

diffusiveFluxBoundaryConditions = {0:{0:getDFBC}}