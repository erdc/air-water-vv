from proteus.ctransportCoefficients import smoothedHeaviside as sH
from proteus.ctransportCoefficients import smoothedHeaviside_integral as sHi
from proteus.AnalyticalSolutions import AS_base

class IC_Base(object):
    def __init__(self, problem=None):
        self.Problem = problem

    @classmethod
    def newGlobalIC(cls, name, default=None):
        """
        Makes a new initial condition with a default value. This creates a new
        class attribute and adds it to all IC_Base class instances.
        :param name: name of the new class attribute
        :param default: default value of new attr (if None: will be a class
                        instance AS_base() of proteus.AnalyticalSolutions)
        """
        setattr(cls, name, default)

class IC_RANS(IC_Base):
    def __init__(self, problem=None):
        super(TwoPhaseFlowIC, self).__init__(problem)
        self.Pressure = AS_base()
        self.VelocityU = AS_base()
        self.VelocityV = AS_base()
        self.VelocityW = AS_base()
        self.VOF = AS_base()
        self.LevelSet = AS_base()
        self.LevelSetConsrv = AS_base()
        self.Redist = AS_base()
        self.Kappa = AS_base()
        self.Dissipation = AS_base()

    def setHorizontalLevelBetweenFluids(self, fluid_level, fluid1='water',
                                        fluid2='air', vert_axis=None,
                                        epsFact_consrv_heaviside=None):
        """
        Sets fluid 0 level in the domain (x is fluid 0 if x[vert_axis] < fluid_level)
        :param fluid_level: level below which fluid is fluid 0
        -----
        This sets/resets the following initial conditions:
        - Pressure
        - LevelSet
        - Redist
        - VOF
        """
        self.initialConditions[fluid1+'Level'] = fluid_level
        pb = self.Problem
        pv = pb.PhysicalVariables
        num = pb.NumericalOptions
        if vert_axis is None:
            vert_axis = pb.nd-1  # y in 2D, z in 3D

        def pressure_init(x, t):
            top = pb.Domain.L[vert_axis]  # 'highest' point in domain
            rho1, rho2, g = fluid1['rho'], fluid2['rho'], pv.g
            p_L = 0.0
            phi_L = top-fluid_level
            phi = x[vert_axis]-fluid_level
            he = pb.Domain.MeshOptions.he
            eps = num.epsFact['consrv_heaviside']
            return p_L-g[vert_axis]*(rho1*(phi_L-phi)+(rho2-rho1) \
                                     *(sHi(eps*he,phi_L)-sHi(eps*he,phi)))
        def phi_init(self, x, t):
            return x[vert_axis]-fluid_level

        def vof_init(self, x, t):
            he = pb.Domain.MeshOptions.he
            eps = num.epsFact['consrv_heaviside']
            return sH(eps*he,x[vert_axis]-fluid_level)

        self.Pressure.uOfXT = pressure_init
        self.LevelSet.uOfXT = phi_init
        self.Redist.uOfXT = phi_init
        self.VOF.uOfXT = vof_init

    def setAtRest(self):
        """
        Sets the initial velocity to 0 for all fluid in the domain
        """
        self.VelocityU.uOfXT = lambda x, t: 0.
        self.VelocityV.uOfXT = lambda x, t: 0.
        self.VelocityW.uOfXT = lambda x, t: 0.
