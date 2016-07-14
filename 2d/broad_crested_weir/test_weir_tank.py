"""
Testing module for weir_tank.py
Based heavily off of test_spatialtools.py
"""
import unittest
import numpy.testing as npt
import numpy as np
from nose.tools import eq_
from proteus import Comm, Profiling
from proteus.Profiling import logEvent
from proteus.Domain import (PiecewiseLinearComplexDomain,
                            PlanarStraightLineGraphDomain)
from proteus.SpatialTools import (Rectangle,
                                  Cuboid,
                                  CustomShape,
                                  assembleDomain)
from proteus.mprans.SpatialTools import (Rectangle as RectangleRANS,
                                         Cuboid as CuboidRANS,
                                         CustomShape as CustomShapeRANS,
                                         assembleDomain as assembleDomainRANS,
                                         Tank2D,
                                         Tank3D,
                                         RigidBody)
from weir_tank import ShapeWeirTank2D

comm = Comm.init()
Profiling.procID = comm.rank()

#[temp] Integration comments for tests/test_spatialtools or similar:
#       Similar names should mostly be filtered into tests of the same name in
#       the classes given there.  There are a few missing (no assembleDomain
#       for instance, because it's only one shape), and a few extras
#       that mostly test minor functions.

logEvent("Testing Weir Tank")

def create_domain2D():
    return PlanarStraightLineGraphDomain()

def create_weir_tank(domain,dim=(0.,0.), obstacle_intersects=None,
                     obstacle_points=None, floating_obstacles=None,
                     floating_centers=None,special_boundaries=None,
                     coords=(0.,0.), from_0=True):
    return ShapeWeirTank2D(domain, dim, obstacle_intersects, obstacle_points,
                           floating_obstacles, floating_centers,
                           special_boundaries, coords, from_0)

def create_elaborate_weir(domain):
    """Builds a complicated geometry"""
    dim = (2.,2.)
    obstacle_intersects = [0,0.9]
    obstacle_points = [[-0.2,-0.9],
                       [0.0,-0.8],
                       [-0.2,-0.7],
                       [-0.1,-0.6],
                       [-0.5,-0.6],
                       [-0.5,-0.5],
                       [0.0,-0.5],
                       [0.0,0.5],
                       [-0.5,0.5],
                       [-0.3,0.6],
                       [0.1,0.5],
                       [0.1,0.0],
                       [0.3,-0.7],
                       [0.05,-0.8],
                       [0.3,-0.9],
                       [0.1,-1.0],
                       [0.4,-1.0],
                       [0.8,-0.9]]
    floating_obstacles = [[[-0.1,-0.4],
                           [-0.8,-0.3],
                           [-0.7,0.3]],
                          [[0.1,0.6],
                           [0.1,0.8],
                           [-0.2,0.8]],
                          [[0.8,-0.1],
                           [0.8,0.2],
                           [0.45,0.2],
                           [0.45,0.45],
                           [0.9,0.45],
                           [0.9,-0.3],
                           [0.45,-0.3],
                           [0.45,-0.1]]]
    special_boundaries = {[0.1,-1.0]: 'any name works'}
    floating_centers = [[-0.5,-0.2],
                        [0.0,0.75],
                        [0.7,-0.2]]
    coords = (0.,0.)
    from_0 = False
    return ShapeWeirTank2D(domain, dim, obstacle_intersects, obstacle_points,
                           floating_obstacles, floating_centers,
                           special_boundaries, coords, from_0)

class TestWeirTank(unittest.TestCase):

    def test_create_shapes(self):
        """
        Test if the shape can be created.
        """
        domain2D = create_domain2D()
        wt = create_weir_tank(domain2D, dim=(2.,2.))
        assert (wt.x0, wt.x1, wt.y0, wt.y1) == (0.,2.,0.,2.)
        wt2 = create_weir_tank(domain2D, dim = (2.,2.), from_0 = False)
        assert (wt2.x0, wt2.x1, wt2.y0, wt2.y1) == (-1., 1., -1., 1.)
        wt3 = create_weir_tank(domain2D, dim = (4.,4.5),
                               from_0 = False, coords = (-110, 110))
        assert (wt3.x0, wt3.x1, wt3.y0, wt3.y1) == (-114, -106, 114.5, 105.5)

    def test_lineVertexInterpolation(self):
        """
        Test if the interpolation method works.
        """
        domain2D = create_domain2D()
        wt = create_weir_tank(domain2D, dim=(2.,2.))
        wt.vertices += [[0,1],[4,2]]
        wt.segments += [[0,1],[1,0]] #both directions tested
        x = 2
        assert 3 == wt._lineVertexInterpolation(x, [0,1])
        assert 3 == wt._lineVertexInterpolation(x, [1,0])
        #[temp] this tests a real basic case.  Make new vertices/segments

    def test_addSegments(self):
        """Test the _addSegments method by adding some segments

        We set up some segments to have special conditions, and some to not
        """
        #[temp]

    def test_pointsOutsideTank(self):
        """Throw a couple points inside and outside the tank - expect error"""
        domain2D = create_domain2D()
        wt = create_weir_tank(domain2D, dim=(1.,1.))
        try:
            wt._checkObstaclePointsWithinTank([0,0])     # corner
            wt._checkObstaclePointsWithinTank([1,1])     # corner
            wt._checkObstaclePointsWithinTank([0.5,0.5]) # middle
            wt._checkObstaclePointsWithinTank([0,0.5])   # edge
        except ValueError:
            self.fail('Obstacle incorrectly declared to be out of tank!')
        self.assertRaises(ValueError,
                          wt._checkObstaclePointsWithinTank, [0,-0.1])   # close
        self.assertRaises(ValueError,
                          wt._checkObstaclePointsWithinTank, [0.5, 100]) # far
        self.assertRaises(ValueError,
                          wt._checkObstaclePointsWithinTank, [100, 0.4]) # far
        self.assertRaises(ValueError,
                          wt._checkObstaclePointsWithinTank, [2, -3])    # both wrong

    def test_constructFloatingObstacle(self):
        """Build a simple floating obstacle and construct."""
        domain2D = create_domain2D()
        wt = create_weir_tank(domain2D,dim=(2.,2.),
                              coords=(0.,0.),from_0=False) #unit square around 0
        wt._constructFloatingObstacleFrame([[-0.5,0],[0.5,-0.5],[0.4,0.3]])
        #[temp] assert there are segments 5-6,6-7,7-8 and that vertices are 5,6,7 equal to the ones above

    def test_initialFrame(self):
        """Test that the initial frame construction works w/ a simple suite"""
        domain2D = create_domain2D()
        wt = create_elaborate_weir(domain2D)
        #[temp]

    def test_setHorizontalVariableMesh(self):
        """Test that the horizontal variable mesh boundaries are set."""
        domain2D = create_domain2D()
        wt = create_weir_tank(domain2D, dim=[2.4,0.5],
                              coords=[0.1,0.], from_0=False) #x0 = -1.1, x1 = 1.3
        wt.setHorizontalVariableMesh([1.2, -1, 0, 0.5, -0.2, 0.6])
        wt.setHorizontalVariableMesh([-1,0.2,-0.5,0.8])
        assert wt.region_boundaries == [-1, -0.5, -0.2, 0, 0.2, 0.5, 0.6, 0.8, 1.2]
        self.assertRaises(ValueError, wt.setHorizontalVariableMesh, [-1.1,])
        self.assertRaises(ValueError, wt.setHorizontalVariableMesh, [1.3, ])
        self.assertRaises(ValueError, wt.setHorizontalVariableMesh, [-1,0.3,4])

    def test_additionalFrame(self):
        """Test that the later frame builds a correct variable region setup"""
        domain2D = create_domain2D()
        wt = create_elaborate_weir(domain2D)
        wt.setHorizontalVariableMesh([-0.5,0.0,-0.1,0.45,0.7])
        #[temp]
