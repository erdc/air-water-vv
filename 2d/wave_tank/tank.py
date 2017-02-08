from proteus import Domain, Context
from proteus.mprans import SpatialTools as st
from proteus import WaveTools as wt
from math import *
import numpy as np



opts=Context.Options([
    # predefined test cases
    ("water_level", 1.5, "Height of free surface above seabed"),
    # tank
    ("tank_dim", (1.596*4, 3.,), "Dimensions of the tank"),
    ("tank_sponge", (1.596*2, 1.596*2), "Length of relaxation zones zones (left, right)"),
    ("tank_BC", 'freeslip', "Length of absorption zones (front/back, left/right)"),
    ("gauge_output", False, "Places Gauges in tank (5 per wavelength)"),
    # waves
    ("waves", True, "Generate waves (True/False)"),
    ("wave_period", 1.003, "Period of the waves"),
    ("wave_height", 0.07, "Height of the waves"),
    ("wave_dir", (1.,0.,0.), "Direction of the waves (from left boundary)"),
    ("wave_wavelength", 1.596, "Direction of the waves (from left boundary)"),
    ("wave_type", 'Fenton', "type of wave"),
    ("Bcoeff", np.array([0.13641547,0.00018183,4.61e-06,1.3e-07,0.0,0.0,0.0,0.0]), "BCoeffs"),
    ("Ycoeff", np.array([0.13676376,0.00961348,0.00102053,0.00012878,1.788e-05,2.64e-06,4.2e-07,1.3e-07]), "YCoeffs"),
    # mesh refinement
    ("refinement", False, "Gradual refinement"),
    ("he", 0.015, "Set characteristic element size"),
    ("he_max", 10, "Set maximum characteristic element size"),
    ("he_max_water", 10, "Set maximum characteristic in water phase"),
    ("refinement_freesurface", 0.1,"Set area of constant refinement around free surface (+/- value)"),
    ("refinement_caisson", 0.,"Set area of constant refinement (Box) around caisson (+/- value)"),
    ("refinement_grading", np.sqrt(1.1*4./np.sqrt(3.))/np.sqrt(1.*4./np.sqrt(3)), "Grading of refinement/coarsening (default: 10% volume)"),
    # numerical options
    ("gen_mesh", True, "True: generate new mesh every time. False: do not generate mesh if file exists"),
    ("use_gmsh", False, "True: use Gmsh. False: use Triangle/Tetgen"),
    ("movingDomain", True, "True/False"),
    ("T", 10.0, "Simulation time"),
    ("dt_init", 0.001, "Initial time step"),
    ("dt_fixed", None, "Fixed (maximum) time step"),
    ("timeIntegration", "backwardEuler", "Time integration scheme (backwardEuler/VBDF)"),
    ("cfl", 0.4 , "Target cfl"),
    ("nsave",  20, "Number of time steps to save per second"),
    ("useRANS", 0, "RANS model"),
    ("parallel", True ,"Run in parallel")])



# ----- CONTEXT ------ #

# waves
if opts.waves is True:
    period = opts.wave_period
    height = opts.wave_height
    mwl = depth = opts.water_level
    direction = opts.wave_dir
    wave = wt.MonochromaticWaves(period=period, waveHeight=height, mwl=mwl, depth=depth,
                                 g=np.array([0., -9.81, 0.]), waveDir=direction,
                                 wavelength=opts.wave_wavelength,
                                 waveType=opts.wave_type,
                                 Ycoeff=np.array(opts.Ycoeff),
                                 Bcoeff=np.array(opts.Bcoeff),
                                 Nf=len(opts.Bcoeff))
    wavelength = wave.wavelength

# tank options
waterLevel = opts.water_level
tank_dim = opts.tank_dim
tank_sponge = opts.tank_sponge

# ----- DOMAIN ----- #

domain = Domain.PlanarStraightLineGraphDomain()
# caisson options


# ----- SHAPES ----- #
tank = st.Tank2D(domain, tank_dim)
tank.setSponge(x_n=tank_sponge[0], x_p=tank_sponge[1])
left = right = False
if tank_sponge[0]: left = True
if tank_sponge[1]: right = True
if left:
    if opts.waves is True:
        smoothing = opts.he*3.
        tank.setGenerationZones(x_n=left, waves=wave, smoothing=smoothing, dragAlpha=0.5/1.004e-6)
        tank.BC['x-'].setUnsteadyTwoPhaseVelocityInlet(wave, smoothing=smoothing, vert_axis=1)
    else:
        tank.setAbsorptionZones(x_n=left, dragAlpha=0.5/1.004e-6)
        tank.BC['x-'].setNoSlip()
if right:
    tank.setAbsorptionZones(x_p=right)
if opts.caisson:
    # let gmsh know that the caisson is IN the tank
    tank.setChildShape(caisson, 0)

# ----- BOUNDARY CONDITIONS ----- #

tank.BC['y+'].setAtmosphere()
if opts.tank_BC == 'noslip':
    tank.BC['y-'].setNoSlip()
if opts.tank_BC == 'freeslip':
    tank.BC['y-'].setFreeSlip()
tank.BC['x+'].setNoSlip()
tank.BC['sponge'].setNonMaterial()

for bc in tank.BC_list:
    bc.setFixedNodes()


# ----- GAUGES ----- #

if opts.gauge_output:
    if left or right:
        gauge_dx = tank_sponge[0]/10.
    else:
        gauge_dx = tank_dim[0]/10.
    probes=np.linspace(-tank_sponge[0], tank_dim[0]+tank_sponge[1], (tank_sponge[0]+tank_dim[0]+tank_sponge[1])/gauge_dx+1)
    PG=[]
    PG2=[]
    LIG = []
    zProbes=waterLevel*0.5
    for i in probes:
        PG.append((i, zProbes, 0.),)
        PG2.append((i, waterLevel, 0.),)
        if i == probes[0]:
            LIG.append(((i, 0.+0.0001, 0.),(i, tank_dim[1]-0.0001,0.)),)
        elif i != probes[0]:
            LIG.append(((i-0.0001, 0.+0.0001, 0.),(i-0.0001, tank_dim[1]-0.0001,0.)),)
    tank.attachPointGauges(
        'twp',
        gauges = ((('p',), PG),),
        activeTime=(0, opts.T),
        sampleRate=0,
        fileName='pointGauge_pressure.csv'
    )
    tank.attachPointGauges(
        'ls',
        gauges = ((('phi',), PG),),
        activeTime=(0, opts.T),
        sampleRate=0,
        fileName='pointGauge_levelset.csv'
    )

    tank.attachLineIntegralGauges(
        'vof',
        gauges=((('vof',), LIG),),
        activeTime = (0., opts.T),
        sampleRate = 0,
        fileName = 'lineGauge.csv'
    )

# ----- ASSEMBLE DOMAIN ----- #

st.assembleDomain(domain)

# ----- REFINEMENT OPTIONS ----- #

if opts.refinement:
    import MeshRefinement as mr
    tank.MeshOptions = mr.MeshOptions(tank)
    grading = opts.refinement_grading
    he2 = opts.he
    def mesh_grading(start, he, grading):
        return '{0}*{2}^(1+log((-1/{2}*(abs({1})-{0})+abs({1}))/{0})/log({2}))'.format(he, start, grading)
    he_max = opts.he_max
    # he_fs = he2
    ecH = 3.
    if opts.refinement_freesurface > 0:
        box = opts.refinement_freesurface
    else:
        box = ecH*he2
    tank.MeshOptions.refineBox(he2, he_max, -tank_sponge[0], tank_dim[0]+tank_sponge[1], waterLevel-box, waterLevel+box)
    tank.MeshOptions.setRefinementFunction(mesh_grading(start='y-{0}'.format(waterLevel-box), he=he2, grading=grading))
    tank.MeshOptions.setRefinementFunction(mesh_grading(start='y-{0}'.format(waterLevel+box), he=he2, grading=grading))
    domain.MeshOptions.LcMax = he_max #coarse grid
    if opts.use_gmsh is True and opts.refinement is True:
        domain.MeshOptions.he = he_max #coarse grid
    else:
        domain.MeshOptions.he = he2 #coarse grid
        domain.MeshOptions.LcMax = he2 #coarse grid
    tank.MeshOptions.refineBox(opts.he_max_water, he_max, -tank_sponge[0], tank_dim[0]+tank_sponge[1], 0., waterLevel)
    mr._assembleRefinementOptions(domain)
    from proteus import Comm
    comm = Comm.get()
    if comm.isMaster():
        mr.writeGeo(domain, 'mesh', append=False)

domain.MeshOptions.use_gmsh = opts.use_gmsh
domain.MeshOptions.genMesh = opts.gen_mesh
domain.use_gmsh = opts.use_gmsh


##########################################
# Numerical Options and other parameters #
##########################################


rho_0=998.2
nu_0 =1.004e-6
rho_1=1.205
nu_1 =1.500e-5
sigma_01=0.0
g = [0., -9.81]




from math import *
from proteus import MeshTools, AuxiliaryVariables
import numpy
import proteus.MeshTools
from proteus import Domain
from proteus.Profiling import logEvent
from proteus.default_n import *
from proteus.ctransportCoefficients import smoothedHeaviside
from proteus.ctransportCoefficients import smoothedHeaviside_integral


#----------------------------------------------------
# Boundary conditions and other flags
#----------------------------------------------------
movingDomain=opts.movingDomain
checkMass=False
applyCorrection=True
applyRedistancing=True
freezeLevelSet=True

#----------------------------------------------------
# Time stepping and velocity
#----------------------------------------------------
weak_bc_penalty_constant = 10.0/nu_0#Re
dt_init = opts.dt_init
T = opts.T
nDTout = int(opts.T*opts.nsave)
timeIntegration = opts.timeIntegration
if nDTout > 0:
    dt_out= (T-dt_init)/nDTout
else:
    dt_out = 0
runCFL = opts.cfl
dt_fixed = opts.dt_fixed

#----------------------------------------------------

#  Discretization -- input options
useOldPETSc=False
useSuperlu = not True
spaceOrder = 1
useHex     = False
useRBLES   = 0.0
useMetrics = 1.0
useVF = 1.0
useOnlyVF = False
useRANS = opts.useRANS # 0 -- None
            # 1 -- K-Epsilon
            # 2 -- K-Omega, 1998
            # 3 -- K-Omega, 1988
# Input checks
if spaceOrder not in [1,2]:
    print "INVALID: spaceOrder" + spaceOrder
    sys.exit()

if useRBLES not in [0.0, 1.0]:
    print "INVALID: useRBLES" + useRBLES
    sys.exit()

if useMetrics not in [0.0, 1.0]:
    print "INVALID: useMetrics"
    sys.exit()

#  Discretization
nd = 2
if spaceOrder == 1:
    hFactor=1.0
    if useHex:
	 basis=C0_AffineLinearOnCubeWithNodalBasis
         elementQuadrature = CubeGaussQuadrature(nd,3)
         elementBoundaryQuadrature = CubeGaussQuadrature(nd-1,3)
    else:
    	 basis=C0_AffineLinearOnSimplexWithNodalBasis
         elementQuadrature = SimplexGaussQuadrature(nd,3)
         elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,3)
         #elementBoundaryQuadrature = SimplexLobattoQuadrature(nd-1,1)
elif spaceOrder == 2:
    hFactor=0.5
    if useHex:
	basis=C0_AffineLagrangeOnCubeWithNodalBasis
        elementQuadrature = CubeGaussQuadrature(nd,4)
        elementBoundaryQuadrature = CubeGaussQuadrature(nd-1,4)
    else:
	basis=C0_AffineQuadraticOnSimplexWithNodalBasis
        elementQuadrature = SimplexGaussQuadrature(nd,4)
        elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,4)


# Numerical parameters
sc = 0.25 # default: 0.5. Test: 0.25
sc_beta = 1. # default: 1.5. Test: 1.
epsFact_consrv_diffusion = 0.1 # default: 1.0. Test: 0.1
ns_forceStrongDirichlet = False
backgroundDiffusionFactor=0.01
if useMetrics:
    ns_shockCapturingFactor  = sc
    ns_lag_shockCapturing = True
    ns_lag_subgridError = True
    ls_shockCapturingFactor  = sc
    ls_lag_shockCapturing = True
    ls_sc_uref  = 1.0
    ls_sc_beta  = sc_beta
    vof_shockCapturingFactor = sc
    vof_lag_shockCapturing = True
    vof_sc_uref = 1.0
    vof_sc_beta = sc_beta
    rd_shockCapturingFactor  =sc
    rd_lag_shockCapturing = False
    epsFact_density    = 3.
    epsFact_viscosity  = epsFact_curvature  = epsFact_vof = epsFact_consrv_heaviside = epsFact_consrv_dirac = epsFact_density
    epsFact_redistance = 0.33
    epsFact_consrv_diffusion = epsFact_consrv_diffusion
    redist_Newton = True#False
    kappa_shockCapturingFactor = sc
    kappa_lag_shockCapturing = False#True
    kappa_sc_uref = 1.0
    kappa_sc_beta = sc_beta
    dissipation_shockCapturingFactor = sc
    dissipation_lag_shockCapturing = False#True
    dissipation_sc_uref = 1.0
    dissipation_sc_beta = sc_beta
else:
    ns_shockCapturingFactor  = 0.9
    ns_lag_shockCapturing = True
    ns_lag_subgridError = True
    ls_shockCapturingFactor  = 0.9
    ls_lag_shockCapturing = True
    ls_sc_uref  = 1.0
    ls_sc_beta  = 1.0
    vof_shockCapturingFactor = 0.9
    vof_lag_shockCapturing = True
    vof_sc_uref  = 1.0
    vof_sc_beta  = 1.0
    rd_shockCapturingFactor  = 0.9
    rd_lag_shockCapturing = False
    epsFact_density    = 1.5
    epsFact_viscosity  = epsFact_curvature  = epsFact_vof = epsFact_consrv_heaviside = epsFact_consrv_dirac = epsFact_density
    epsFact_redistance = 0.33
    epsFact_consrv_diffusion = 10.0
    redist_Newton = False#True
    kappa_shockCapturingFactor = 0.9
    kappa_lag_shockCapturing = True#False
    kappa_sc_uref  = 1.0
    kappa_sc_beta  = 1.0
    dissipation_shockCapturingFactor = 0.9
    dissipation_lag_shockCapturing = True#False
    dissipation_sc_uref  = 1.0
    dissipation_sc_beta  = 1.0

ns_nl_atol_res = 1e-6 #max(1.0e-12,tolfac*he**2)
vof_nl_atol_res = 1e-6 #max(1.0e-12,tolfac*he**2)
ls_nl_atol_res = 1e-6 #max(1.0e-12,tolfac*he**2)
mcorr_nl_atol_res = 1e-6 #max(1.0e-12,0.1*tolfac*he**2)
rd_nl_atol_res = 1e-4 #max(1.0e-12,tolfac*he)
kappa_nl_atol_res = 1e-6 #max(1.0e-12,tolfac*he**2)
dissipation_nl_atol_res = 1e-6 #max(1.0e-12,tolfac*he**2)
mesh_nl_atol_res = 1e-6 #max(1.0e-12,opts.mesh_tol*he**2)

#turbulence
ns_closure=0 #1-classic smagorinsky, 2-dynamic smagorinsky, 3 -- k-epsilon, 4 -- k-omega

if useRANS == 1:
    ns_closure = 3
elif useRANS >= 2:
    ns_closure == 4

def twpflowPressure_init(x, t):
    p_L = 0.0
    phi_L = tank_dim[nd-1] - waterLevel
    phi = x[nd-1] - waterLevel
    return p_L -g[nd-1]*(rho_0*(phi_L - phi)+(rho_1 -rho_0)*(smoothedHeaviside_integral(epsFact_consrv_heaviside*opts.he,phi_L)
                                                         -smoothedHeaviside_integral(epsFact_consrv_heaviside*opts.he,phi)))
