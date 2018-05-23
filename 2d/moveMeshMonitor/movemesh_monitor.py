import matplotlib.pyplot as plt
import numpy as np
from proteus.iproteus import *
import proteus.default_p as physics
import proteus.default_n as numerics
from proteus.TransportCoefficients import PoissonEquationCoefficients
from proteus  import Profiling
from proteus.Profiling import logEvent
Profiling.logLevel=7
Profiling.verbose=True


def get_center_area(e_nodes):
    detJ = (e_nodes[1][0] - e_nodes[0][0]) * (e_nodes[2][1] -
                                              e_nodes[0][1]) - (e_nodes[1][1] - e_nodes[0][1]) * (e_nodes[2][0] -
                                                                                                  e_nodes[0][0])
    # since the orientation is clockwise
    center = ( e_nodes[0]+e_nodes[1]+e_nodes[2] )/3.
    area = 0.5 * np.abs(detJ)
    return area, center

import copy

he_max=1.
he_min=0.025
r = 0.1
nSmooth = 10



nd = 2

# use structured mesh
# domain =  Domain.RectangularDomain(L=[1.0,1.0], x=[0.,0.])
# nn = 50
from proteus import SpatialTools as st
domain = Domain.PlanarStraightLineGraphDomain()
rect = st.Rectangle(domain, dim=[1.,1.])
rect.translate(np.array([rect.dim[0]/2., rect.dim[1]/2.]))
boundaryNormals = {rect.boundaryTags['x-']: np.array([-1.,0.]),
                   rect.boundaryTags['x+']: np.array([1.,0.]),
                   rect.boundaryTags['y-']: np.array([0.,-1.]),
                   rect.boundaryTags['y+']: np.array([0.,1.])}
boundaryNormals = None

# center1 = [1.5-2*r, 1.5-0.1]
# center1 = [rect.dim[0]/2., rect.dim[1]/2.+0.2]
# center2 = [rect.dim[0]/2., rect.dim[1]/2.-0.2]
# center3 = [rect.dim[0]/2.+0.2, rect.dim[1]/2.]
# center4 = [rect.dim[0]/2.-0.2, rect.dim[1]/2.]

center1 = [rect.dim[0]/2., rect.dim[1]/2.]

st.assembleDomain(domain)
he = 1./50.
domain.MeshOptions.he = he
domain.MeshOptions.setTriangleOptions()
domain.MeshOptions.setOutputFiles(name="mesh")
domain.use_gmsh = False
genMesh = True

use_gmsh = False
if use_gmsh:
    from MeshRefinement import geometry_to_gmsh
    import py2gmsh
    domain.use_gmsh = True
    mesh = geometry_to_gmsh(domain)
    field_list = []
    me01 = py2gmsh.Fields.MathEval(mesh=mesh)
    me01.F = "{he}+{he}*10*abs(sqrt((x-{center_x})^2+(y-{center_y})^2)-{radius})".format(he=he, center_x=center1[0], center_y=center1[1], radius=r)
    field_list += [me01]
    # me02 = py2gmsh.Fields.MathEval(mesh=mesh)
    # me02.F = "{he}+{he}*10*abs(sqrt((x-{center_x})^2+(y-{center_y})^2)-{radius})".format(he=he, center_x=center2[0], center_y=center2[1], radius=r)
    # field_list += [me02]
    # me03 = py2gmsh.Fields.MathEval(mesh=mesh)
    # me03.F = "{he}+{he}*10*abs(sqrt((x-{center_x})^2+(y-{center_y})^2)-{radius})".format(he=he, center_x=center3[0], center_y=center3[1], radius=r)
    # field_list += [me03]
    # me04 = py2gmsh.Fields.MathEval(mesh=mesh)
    # me04.F = "{he}+{he}*10*abs(sqrt((x-{center_x})^2+(y-{center_y})^2)-{radius})".format(he=he, center_x=center4[0], center_y=center4[1], radius=r)
    # field_list += [me04]
    fmin = py2gmsh.Fields.Min(mesh=mesh)
    fmin.FieldsList = field_list
    mesh.setBackgroundField(fmin)
    mesh.Options.Mesh.CharacteristicLengthMax = he*(he_max/he_min)
    mesh.Options.Mesh.CharacteristicLengthMax = he*5
    mesh.writeGeo("mesh.geo")
    domain.geofile = "mesh"


def my_func(x):
    return min(he_max, max(np.abs(np.sqrt((x[0]-0.5)**2+(x[1]-0.5)**2)-r), he_min))
scale = 0.5
def my_func(x, t):
    t_start = 1000.
    if t>t_start:
        dist1 = np.sqrt((x[0]-(center1[0]))**2+(x[1]-(center1[1]-(t-t_start)/100))**2)-r
        # dist2 = np.sqrt((x[0]-(center2[0]))**2+(x[1]-(center2[1]+(t-t_start)/100))**2)-r
        # dist3 = np.sqrt((x[0]-(center3[0])+(t-t_start)/100)**2+(x[1]-(center3[1]))**2)-r
        # dist4 = np.sqrt((x[0]-(center4[0])-(t-t_start)/100)**2+(x[1]-(center4[1]))**2)-r
        circle1 = abs(dist1)/scale
        # circle2 = min(he_max, max(abs(dist2)/scale, he_min))
        # circle3 = min(he_max, max(abs(dist3)/scale, he_min))
        # circle4 = min(he_max, max(abs(dist4)/scale, he_min))
    else:
        dist1 = np.sqrt((x[0]-center1[0])**2+(x[1]-center1[1])**2)-r
        # dist2 = np.sqrt((x[0]-center2[0])**2+(x[1]-center2[1])**2)-r
        # dist3 = np.sqrt((x[0]-center3[0])**2+(x[1]-center3[1])**2)-r
        # dist4 = np.sqrt((x[0]-center4[0])**2+(x[1]-center4[1])**2)-r
        circle1 = abs(dist1)/scale
        # circle2 = min(he_max, max(abs(dist2)/scale, he_min))
        # circle3 = min(he_max, max(abs(dist3)/scale, he_min))
        # circle4 = min(he_max, max(abs(dist4)/scale, he_min))
    #border = min(he_max, min(max(abs(x[0]-0.)/scale, he_min),
    #                         max(abs(x[0]-domain.L[0])/scale, he_min),
    #                         max(abs(x[1]-0.)/scale, he_min),
    #                         max(abs(x[1]-domain.L[1])/scale, he_min)))
    # return min(circle1, circle2, circle3, circle4) #min(border, circle)
    return circle1#min(circle1, circle2) #min(border, circle)


movingDomain = True




from petsc4py import PETSc
OptDB = PETSc.Options()
#OptDB.setValue("ksp_type", "cg")
#OptDB.setValue("pc_type", "jacobi")
#OptDB.setValue("pc_factor_mat_solver_package", "superlu_dist")
OptDB.setValue("ksp_constant_null_space", 1)
#OptDB.setValue("ksp_converged_true_residual", 1)
OptDB.setValue("info", 1)
#OptDB.setValue("pc_factor_shift_type","NONZERO")
#OptDB.setValue("pc_factor_shift_amount", 1e-10)
print(OptDB.getAll())



from proteus.default_n import *
movingDomain = True
T = 10.
tnList = [i*0.1 for i in range(1000)]
spaceOrder = 1
useHex = False
nd = domain.nd
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


# model = ns.modelList[0].levelModelList[-1]

# x0 = model.mesh_array0[:,0]
# y0 = model.mesh_array0[:,1]
# x1 = model.mesh.nodeArray[:,0]
# y1 = model.mesh.nodeArray[:,1]
# z = model.u[0].dof
# qx0 = model.q['x'][:,:,0].flatten()
# qy0 = model.q['x'][:,:,1].flatten()
# qzx = model.q[('grad(u)',0)][:,:,0].flatten()
# qzy = model.q[('grad(u)',0)][:,:,1].flatten()
# Nlevels = 1000

# import matplotlib.cm as cm

# f, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(10,5))
# CS = ax1.tricontourf(x0,
#                 y0,
#                 model.mesh.elementNodesArray,
#                 z,
#                 Nlevels,
#                 cmap=cm.jet)
# ax1.set_aspect('equal')
# ax2.tricontourf(x1,
#                 y1,
#                 model.mesh.elementNodesArray,
#                 z,
#                 Nlevels,
#                 cmap=cm.jet)
# ax2.set_aspect('equal')
# print('min/max',min(model.u[0].dof), max(model.u[0].dof))
# c0 = CS.collections[0]
# f.colorbar(CS, ax=ax1)
# f.savefig('sol.png', bbox_inches='tight', dpi=100)

# # plt.figure()
# # plt.plot(model.q['x'][:,:,0], model.q['x'][:,:,1], 'b.')

# import matplotlib.tri as tri
# # Create the Triangulation; no triangles so Delaunay triangulation created.
# f, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(10,5))
# ax1.triplot(x0, y0, model.mesh.elementNodesArray, lw=0.8)
# ax1.set_aspect('equal')
# ax2.triplot(x1, y1, model.mesh.elementNodesArray, lw=0.8)
# ax2.set_aspect('equal')
# f.savefig('mesh.png', bbox_inches='tight', dpi=100)

# from mpl_toolkits.mplot3d import Axes3D
# f, (ax1, ax2) = plt.subplots(1, 2, sharey=True,subplot_kw=dict(projection='3d'))
# fig = plt.figure()
# #ax = fig.add_subplot(111, projection='3d')
# ax1.plot_trisurf(x0, y0, z)
# ax2.plot_trisurf(x1, y1, z)
# f, (ax1) = plt.subplots(1, 1, sharey=True, figsize=(10,10), frameon=False)
# CS = ax1.tricontourf(x0,
#                 y0,
#                 model.mesh.elementNodesArray,
#                 z,
#                 Nlevels,
#                 cmap=cm.jet)
# ax1.set_axis_off()
# ax1.set_aspect('equal')
# f.savefig('a.png', bbox_inches='tight', dpi=100)

# f, (ax1) = plt.subplots(1, 1, sharey=True, figsize=(10,10), frameon=False)
# CS = ax1.tricontourf(x0,
#                 y0,
#                 model.mesh.elementNodesArray,
#                 model.grads[:,0],
#                 cmap=cm.jet)
# f.savefig('gradx.png', bbox_inches='tight', dpi=100)

# f, (ax1) = plt.subplots(1, 1, sharey=True, figsize=(10,10), frameon=False)
# CS = ax1.tricontourf(x0,
#                     y0,
#                     model.mesh.elementNodesArray,
#                     model.grads[:,1],
#                     Nlevels,
#                     cmap=cm.jet)
# f.savefig('grady.png', bbox_inches='tight', dpi=100)


# def my_func(x, y):
#     return np.minimum(he_max, np.maximum(np.abs(np.sqrt((x-0.5)**2+(y-0.5)**2)-0.25)/0.25, he_min))
# x = np.linspace(0., 1., 100)
# y = np.linspace(0., 1., 100)
# X, Y = np.meshgrid(x, y)
# Z = my_func(X, Y)
# f, (ax1) = plt.subplots(1, 1, sharey=True, figsize=(10,10))
# CS = ax1.contourf(X,
#                     Y,
#                     Z,
#                     Nlevels,
#                     cmap=cm.jet)
# f.colorbar(CS, ax=ax1)
# f.savefig('func.png', bbox_inches='tight', dpi=100)
