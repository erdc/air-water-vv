
from proteus.default_p import *
from proteus.mprans import MoveMesh
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


nMediaTypes = len(domain.regionFlags)  # (!) should be region flags
smTypes     = np.zeros((nMediaTypes+1, 2), 'd')
smFlags     = np.zeros((nMediaTypes+1,), 'i')
smTypes[:, 0] = 1.0    # E
smTypes[:, 1] = 0.3    # nu

el_max = 1.
el_min = 0.1

def my_func2(x):
    return 2*x[1]*(1-x[1])+2*x[0]*(1-x[0])
class AnalyticalSolution(object):
    def __init__(self):
        pass
    def uOfXT(self, x, t):
        return x[1]*(1-x[1])*x[0]*(1-x[0])
class MyFunc(object):
    def __init__(self):
        self.C = 1.
    def my_func(self, x):
        return min(el_max, max(np.abs(np.sqrt((x[0]-0.5)**2+(x[1]-0.5)**2)-0.25)/0.25, el_min))
    def my_func_scaled(self, x):
        return self.my_func(x)*self.C
    def my_func_reciprocal(self, x):
        return 1./self.my_func(x)
    def my_func_reciprocal_scaled(self, x):
        return 1./self.my_func_scaled(x)
    def RHS(self, x):
        return my_func2(x)
        return 1.
        return self.my_func_reciprocal_scaled(x)-1


myfuncclass = MyFunc()

AS = AnalyticalSolution()
# no LevelModel type
analyticalSolution = {0: AS}
aOfX = {0: lambda x: np.array([[1,0], [0,1]])}
fOfX = {0: my_func2}
nc = 1
nd = domain.nd


def get_center_area(e_nodes):
    detJ = (e_nodes[1][0] - e_nodes[0][0]) * (e_nodes[2][1] -
                                              e_nodes[0][1]) - (e_nodes[1][1] - e_nodes[0][1]) * (e_nodes[2][0] -
                                                                                                  e_nodes[0][0])
    # since the orientation is clockwise
    center = ( e_nodes[0]+e_nodes[1]+e_nodes[2] )/3.
    area = 0.5 * np.abs(detJ)
    return area, center

from proteus import TransportCoefficients
class MyCoeff(TransportCoefficients.PoissonEquationCoefficients):
    def __init__(self):
        # super(MyCoeff, self).__init__(aOfXT, fOfX, nc, nd)
        TransportCoefficients.PoissonEquationCoefficients.__init__(self, aOfX, fOfX, nc, nd)

    def preStep(self, t, firstStep=False):
        integral_1_over_f = 0.0
        
        for e in xrange(self.mesh.elementNodesArray.shape[0]):
            area, center = get_center_area(self.mesh.nodeArray[self.mesh.elementNodesArray[e]])
            integral_1_over_f += area/myfuncclass.my_func(center)
        myfuncclass.C = integral_1_over_f

#        for node in self.mesh.nodeArray:
#            if node[0] < 1e-5 and node[1] < 1e-6:
#                print(node)
#        

    def attachModels(self,modelList):
        """
        Give the TC object access to other models in a loosely coupled split operator formulation (e.g. a transport equation for concentration might get velocity from a flow equation)
        """
        self.model = modelList[-1]

    def initializeMesh(self, mesh):
        self.mesh = mesh


#coefficients = TransportCoefficients.PoissonEquationCoefficients(aOfX,fOfX,nc,nd)
coefficients = MyCoeff()

dirichletConditions = {0: lambda x, flag: domain.bc[flag].hx_dirichlet.init_cython()}

fluxBoundaryConditions = {0: 'noFlow',
                          1: 'noFlow'}

advectiveFluxBoundaryConditions = {}

def getDiffFluxBC5(x,flag):
    if flag != 0: 
        n = numpy.zeros((nd,),'d'); n[0]=1.0
        return lambda x,t: 0.

diffusiveFluxBoundaryConditions = {0: {0: getDiffFluxBC5}}

if nd == 3:
    dirichletConditions[2] = lambda x, flag: domain.bc[flag].hz_dirichlet.init_cython()
    fluxBoundaryConditions[2] = 'noFlow'
    diffusiveFluxBoundaryConditions[2] = {}
    stressFluxBoundaryConditions[3] = lambda x, flag: domain.bc[flag].w_stress.init_cython()

#########################
class My_OneLevelTransport(OneLevelTransport):
    def getJacobian(self,jacobian,skipMassTerms=False):
	print "I am here to assemble Jacobian"
	OneLevelTransport.getJacobian(self,jacobian,skipMassTerms)
	rowptr, colind, nzval = jacobian.getCSRrepresentation()
	for j in range(rowptr[0],rowptr[1]):
		if 0 == colind[j]:
			nzval[j] = 10e20
		else:
			nzval[j] = 0.0	

	


LeveoModelType = My_OneLevelTransport

