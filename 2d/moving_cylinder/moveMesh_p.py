from proteus.default_p import *
from proteus.mprans import MoveMesh
from proteus import Context, AuxiliaryVariables
ct = Context.get()
genMesh = ct.genMesh
movingDomain = ct.movingDomain
L = ct.L
T = ct.T
nd = ct.nd
domain = ct.domain


initialConditions = None

analyticalSolution = {}

nMediaTypes=1
smTypes      = numpy.zeros((nMediaTypes+1,2),'d')
smFlags      = numpy.zeros((nMediaTypes+1,),'i')

smTypes[0,0] = 1.0    ##E
smTypes[0,1] = 0.3    ##nu
smTypes[1,0] = 1.0    ##E
smTypes[1,1] = 0.3    ##nu

LevelModelType = MoveMesh.LevelModel
coefficients = MoveMesh.Coefficients(nd=ct.nd,
				     V_model=1,
                                     modelType_block=smFlags,
				     modelParams_block=smTypes,
                                     meIndex=0)

class FloatingObstacle(AuxiliaryVariables.AV_base):
    import numpy as np
    def __init__(self):
        self.object = None
    def attachModel(self,model,ar):
        self.model=model
        return self
    def attachAuxiliaryVariables(self,avDict):
        self.object = avDict['twp_navier_stokes_p'][0]
    def hx(self,x,t):
        if self.object == None:
            return 0.0
        else:
            #hx = self.object.body.getRelPointPos(np.dot(self.object.last_rotation_inv,(x-self.object.last_position)))[0] - x[0]
            hx = self.object.h[0]
            return hx
    def hy(self,x,t):
        if self.object == None:
            return 0.0
        else:
            #hy = self.object.body.getRelPointPos(np.dot(self.object.last_rotation_inv,(x-self.object.last_position)))[1] - x[1]
            hy = self.object.h[1]
            return hy
    def calculate(self):
        pass

fo = FloatingObstacle()

def getDBC_hx(x,flag):
    if flag in [ct.boundaryTags['left'],
                ct.boundaryTags['right'],
                ct.boundaryTags['top'],
                ct.boundaryTags['bottom'],
                ct.boundaryTags['front'],
                ct.boundaryTags['back']]:
        return lambda x,t: fo.hx(x,t)
    if flag == ct.boundaryTags['obstacle']:
        return lambda x,t: fo.hx(x,t)

def getDBC_hy(x,flag):
    if flag in [ct.boundaryTags['left'],
                ct.boundaryTags['right'],
                ct.boundaryTags['top'],
                ct.boundaryTags['bottom'],
                ct.boundaryTags['front'],
                ct.boundaryTags['back']]:
        return lambda x,t: fo.hy(x,t)
    if flag == ct.boundaryTags['obstacle']:
        return lambda x,t: fo.hy(x,t)

dirichletConditions = {0:getDBC_hx,
                       1:getDBC_hy}

advectiveFluxBoundaryConditions =  {}

diffusiveFluxBoundaryConditions = {0:{},
                                   1:{}}

def stress_u(x,flag):
    return None

def stress_v(x,flag):
    return None

stressFluxBoundaryConditions = {0:stress_u,
                                1:stress_v}
