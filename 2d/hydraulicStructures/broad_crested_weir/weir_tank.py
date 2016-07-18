"""Tools to build a Weir-Tank model in Proteus.

This module sets up the geometric shape and physical regions of interest
of a Broad or Sharp Crested Weir problem, using MPRANs SpatialTools and
BoundaryConditions modules, but may apply to other similar problems.

A 2d weir-tank problem involves a rectangular tank modeled in
two dimensions containing some mixture of air and water which flows
through it and interacts with some form of static obstacle or obstacles.

These obstacles must be representable as some conglomeration of 
polygonal elements, one of which may be attached directly to the floor
of the tank (note that this can be used to produce multiple objects connected
to the tank, or to create a tank bottom which is not flat).

It may include generation and absorption zones outside the tank, horizontally
varying refinement levels, the wave motions as formulated in the WaveTools
module, and the inclusion of special boundary zones on the obstacle
(used, for instance, to place air vents on a weir).
"""

# [temp] Tests Needed:
# ----------------------------------
# Creation of tank
# Set dimensions
# Segment/vertex logic
# Region logic
# Hole logic
# Adding in fake segments
# Setting up horizontal variable mesh (that the variables stick)
# Variable Mesh used.
# Complex geometry tests???

import numpy as np
from proteus import AuxiliaryVariables, Comm, Gauges
from proteus.Profiling import logEvent
from proteus.mprans import BoundaryConditions as BC
from proteus.mprans.SpatialTools import Tank2D
#[temp] inheriting from ShapeRANS is another, potentially more stable, option
#[temp] - but it could require additional work due to needing to work around
#[temp] the generic dimensions.


# [temp] refactor imports at the end to avoid having too many


class TankWithObstacles2D(Tank2D):
    """Build the Weir-Tank geometry & regions of interest.
    
    Parameters
    ----------
    domain: proteus.Domain.D_base
        Domain class instance that holds all the geometrical information
        and boundary conditions of the shape.
    dim: Optional[array_like]
        Dimensions of the rectangular tank observed. (x,y)
    obstacle_intersects: Optional[array_like]
        The x components of two coordinates where the obstacle, or
        weir, connects to the -y edge of the tank.
    obstacle_points: Optional[array_like]
        An array of (x,y) coordinates outlining the shape of the obstacle.
        The obstacle is traced out in a left-to-right manner.  That is, the
        class starts from the smaller obstacle_intersects coordinate and draws
        a segment connecting to the first element of obstacle_points, then
        draws segments connecting each point in the order given, until drawing
        a final segment to the larger obstacle_intersects coordinate.
    floating_obstacles: Optional[array_like]
        An array of arrays of (x,y) coordinates outlining the shapes of any
        obstacle parts that are not connected to the floor.  Each sub array is
        traced out in a left-to-right manner as with obstacle points, except
        that it starts and ends on the first element (closure).
        These obstacles must not overlap.
    floating_centers: Optional[array_like]
        An array of (x,y) coordinates of a point inside each floating obstacle
        previously defined, in order to determine which parts of the region
        to empty.
    special_boundaries: Optional[dict]
        A dictionary of special boundary conditions by name, keyed to a
        2-element list of vertices (by coordinates, not by index) which match
        the start and end of the segment (not necessarily in order).
    coords: Optional[array_like]
        Coordinates of the centroid of the tank (ignoring any obstacles).
        Defaults to 0,0.  Overriden by from_0 option.
    from_0: Optional[boolean]
        If True (default), the tank extends from the origin to positive x,y,z.
        If False, then the tank is centered around coords.
    """  # [deptemp] a later extension which would change fake-segment calculation would be to change the bottom of the tank via some API I can't think up of yet. It would just need to place the intersections intelligently (and the sponge layers might be a problem too...) based on the slopes of the bottom.  Later tests will need this, but not yet, and it's a clear extension off of this.
    count = 0

    def __init__(self, domain, dim=(0., 0.),
                 obstacle_intersects=None, obstacle_points=None,
                 floating_obstacles=None, floating_centers=None,
                 special_boundaries=[],
                 coords=(0., 0.), from_0=True):
        super(TankWithObstacles2D, self).__init__(domain)
        self.__class__.count += 1
        self.name = '2D_weir' + str(self.__class__.count)
        self.dim = L, H = np.array(dim)
        self.obstacle_intersects = obstacle_intersects
        self.obstacle_points = obstacle_points
        self.floating_obstacles = floating_obstacles
        self.special_boundaries = special_boundaries
        self.floating_centers = floating_centers
        self.from_0 = from_0
        self.gauges = {}

        self.spongeLayers = {'x-': None,
                             'x+': None}

        self.region_column = {}
        self.region_boundaries = []

        self.dim = L, H = dim
        if self.from_0 is True:
            x, y = L / 2., H / 2.
        else:
            x, y = coords
        self.coords = [x, y]
        self.x0, self.x1 = x - 0.5 * L, x + 0.5 * L
        self.y0, self.y1 = y - 0.5 * H, y + 0.5 * H

        self.boundaryTags = {'y-': 1, 'x+': 2, 'y+': 3, 'x-': 4,
                             'sponge': 5, 'empty': 6}
        self.b_or = np.array([[0., -1.],
                              [1., 0.],
                              [0., 1.],
                              [-1., 0.]])
        self.BC = {'y-': self.BC_class(shape=self, name='y-',
                                       b_or=self.b_or, b_i=0),
                   'x+': self.BC_class(shape=self, name='x+',
                                       b_or=self.b_or, b_i=1),
                   'y+': self.BC_class(shape=self, name='y+',
                                       b_or=self.b_or, b_i=2),
                   'x-': self.BC_class(shape=self, name='x-',
                                       b_or=self.b_or, b_i=3),
                   'sponge': self.BC_class(shape=self, name='sponge'),
                   'empty': self.BC_class(shape=self, name='empty')}

        next_tag_index = len(self.BC)
        for x in self.special_boundaries:
            if self.special_boundaries[x] not in self.BC:
                name = self.special_boundaries[x]
                self.BC[name] = self.BC_class(shape=self,name=name)
                self.boundaryTags[name] = next_tag_index
                next_tag_index += 1

        self.BC_list = self.BC.values()
        for i in range(4):
             self.BC_list[i].setTank()

        self.BC['y-'].setFreeSlip()
        self.BC['y+'].setAtmosphere()

        self.BC['x-'].setFreeSlip()
        self.BC['x+'].setFreeSlip()

        if obstacle_intersects and not obstacle_points:
            raise ValueError('ERROR: Obstacle intersects defined,'
                             + ' but no points are given.')
        elif obstacle_points and not obstacle_intersects:
            raise ValueError('ERROR: Floor-connected obstacle defined, '
                             + 'but no intercept points given.')
        if floating_obstacles and floating_centers:
            if len(floating_obstacles) != len(floating_centers):
                raise ValueError('ERROR: Number of floating obstacles is'
                                 + ' inconsistent between inputs.')
        elif floating_obstacles and not floating_centers:
            raise ValueError('ERROR: Floating obstacles are not given'
                             + ' hole locations. Cannot parse mesh.')
        elif floating_centers and not floating_obstacles:
            raise ValueError('ERROR: Holes are given without obstacles.')


    def _lineVertexInterpolation(self, x, segment):
        """Given x, returns corresponding y on the line.

        Parameters
        ----------
        x: number
            x coordinate of point desired
        segment: array_like
            segment describing a line you want to put the point on
            (two vertex indices)

        Returns
        -------
        y: number
            corresponding y coordinate (to x) of point on the line
        """
        segment_vertices = [self.vertices[segment[0]],
                            self.vertices[segment[1]]]
        segment_vertices.sort()
        # A very naive algorithm - could cause numerical issues later on.
        run = segment_vertices[1][0] - segment_vertices[0][0]
        rise = segment_vertices[1][1] - segment_vertices[0][1]
        x_run = x - segment_vertices[0][0]
        x_rise = x_run * rise / run
        return segment_vertices[0][1] + x_rise

    def _addSegments(self, segments, segmentFlags, override = False):
        """Add in segments while checking for special boundaries.

        Checks against user defined special boundaries while attaching a list
        of segments to the shape.  It will follow the associated segmentFlag
        if no special boundary condition is given, but will replace this
        with a user defined special flag if this option was selected.

        It will also propagate the change of segment boundary to both endpoints
        of the segment, but not to any other vertices which may lie on it.

        An optional override is included to allow the segments to be added in a
        consistent manner (and one which will update if further conditions or
        tests are included in this method) without the additional cost and
        potential mistakes of special boundary checking (if you know a segment
        has a clear boundary - such as splitting up segments after special
        boundaries have already been placed)

        Parameters
        ----------
        segments: array_like
        segmentFlags: array_like
        override: Optional[boolean]
            If set to True (False is default) this
        """
        if self.segments is None:
            self.segments = []
        if self.segmentFlags is None:
            self.segmentFlags = []
        if not override:
            for i in range(len(segments)):
                self.segments += [segments[i],]
                desired_flag = segmentFlags[i]
                index_0, index_1 = segments[i]
                test1 = (self.vertices[index_0], self.vertices[index_1])
                test2 = (self.vertices[index_1], self.vertices[index_0])
                if test1 in self.special_boundaries:
                    desired_flag = self.boundaryTags[self.special_boundaries[test1]]
                elif test2 in self.special_boundaries:
                    desired_flag = self.boundaryTags[self.special_boundaries[test2]]
                self.vertexFlags[index_0] = desired_flag
                self.vertexFlags[index_1] = desired_flag
                self.segmentFlags += [desired_flag,]
        else:
            self.segments += segments
            self.segmentFlags += segmentFlags


    def _checkObstaclePointsWithinTank(self, vertex):
        """Given an x,y coordinate, check if it is inside the tank.

        Parameters
        ----------
        vertex: array_like
        """
        if not self.x0 <= vertex[0] <= self.x1 or \
                not self.y0 <= vertex[1] <= self.y1:
            raise ValueError(
                'ERROR: ' + str(vertex) + ' is out of tank bounds: '
                + '[{0},{1}] to [{2},{3]]'.format(str(self.x0), str(self.y0),
                                                  str(self.x1), str(self.y1))
            )

    def _constructFloatingObstacleFrame(self, vertex_list):
        """Set up the segment frame for one floating obstacle.

        Parameters
        ----------
        vertex_list: array_like
            A list of (x,y) coordinates to construct the obstacle
            with.  The obstacle is created by tracing a segment
            between each two vertices in order, then wrapping back
            to the first vector in the end.
        """
        new_segments = []
        new_segmentFlags = []
        for vertex in vertex_list:
            self._checkObstaclePointsWithinTank(vertex)
            self.vertices += [vertex, ]
            self.vertexFlags += self.boundaryTags['y-']
            if vertex_list.index(vertex) != len(vertex_list) - 1:
                new_segments += [len(self.vertices) - 1,
                                 len(self.vertices)]
                new_segmentFlags += self.boundaryTags['y-']
            else:
                new_segments += [len(self.vertices) - 1,
                                 len(self.vertices) - len(vertex_list)]
                new_segmentFlags += self.boundaryTags['y-']
        self._addSegments(new_segments, new_segmentFlags)

    def _constructInitialFrame(self):
        """Mark the physical boundaries of the tank-weir problem

        The segments and vertices that form physical geometry for
        the tank are set. 
        """
        leftSponge = self.spongeLayers['x-'] or 0
        rightSponge = self.spongeLayers['x+'] or 0
        obstacle_points = self.obstacle_points or []

        self.vertices = [[self.x0, self.y0],  # 0
                         [self.x1, self.y0],  # 1
                         [self.x1, self.y1],  # 2
                         [self.x0, self.y1]]  # 3
        self.vertexFlags = [self.boundaryTags['y-'],
                            self.boundaryTags['y+'],
                            self.boundaryTags['y+'],
                            self.boundaryTags['y-']]
        if self.obstacle_intersects:
            self.obstacle_intersects.sort()
            self.vertices += [[self.obstacle_intersects[0], self.y0],
                              [self.obstacle_intersects[1], self.y0]]
            self.vertexFlags += [self.boundaryTags['y-'],
                                 self.boundaryTags['y-']]
            new_segments = [[0, 4], [5, 1]]
            new_segmentFlags = [self.boundaryTags['y-'],
                                self.boundaryTags['y-']]
            start_vertex = 4
            end_vertex = 5
            for point in obstacle_points:
                self._checkObstaclePointsWithinTank(point)
                self.vertices += [point, ]
                self.vertexFlags += [self.boundaryTags['y-'], ]
                new_vertex = len(self.vertices) - 1
                new_segments += [[start_vertex, new_vertex], ]
                new_segmentFlags += [self.boundaryTags['y-'], ]
                start_vertex = new_vertex
            new_segments += [[start_vertex, end_vertex], ]
            new_segmentFlags += [self.boundaryTags['y-'], ]
        else:
            new_segments = [[0, 1], ]
            new_segmentFlags = [self.boundaryTags['y-'], ]
        if leftSponge:
            self.vertices += [[self.x0 - leftSponge, self.y0],
                              [self.x0 - leftSponge, self.y1]]
            self.vertexFlags += [self.boundaryTags['y-'],
                                 self.boundaryTags['y+']]
            ls_y0 = len(self.vertices) - 2
            ls_y1 = len(self.vertices) - 1
            new_segments += [[ls_y0, 0],
                             [ls_y1, 3],
                             [ls_y0, ls_y1],
                             [0, 3]]
            new_segmentFlags += [self.boundaryTags['y-'],
                                 self.boundaryTags['y+'],
                                 self.boundaryTags['x-'],
                                 self.boundaryTags['sponge']]
        else:
            new_segments += [[0, 3], ]
            new_segmentFlags += [self.boundaryTags['x-'], ]
        if rightSponge:
            self.vertices += [[self.x1 + rightSponge, self.y0],
                              [self.x1 + rightSponge, self.y1]]
            self.vertexFlags += [self.boundaryTags['y-'],
                                 self.boundaryTags['y+']]
            rs_y0 = len(self.vertices) - 2
            rs_y1 = len(self.vertices) - 1
            new_segments += [[1, rs_y0],
                             [2, rs_y1],
                             [rs_y0, rs_y1],
                             [1, 2]]
            new_segmentFlags += [self.boundaryTags['y-'],
                                 self.boundaryTags['y+'],
                                 self.boundaryTags['x+'],
                                 self.boundaryTags['sponge']]
        else:
            new_segments += [[1, 2], ]
            new_segmentFlags += [self.boundaryTags['x+'], ]
        if self.floating_obstacles:
            for obstacle in self.floating_obstacles:
                self._constructFloatingObstacleFrame(obstacle)
        self._addSegments(new_segments, new_segmentFlags)

    def _constructAdditionalFramework(self):
        """Refine the boundaries of the tank to fit user defined regions

        Construct a new set of vertices and segments that divide
        the physical segmenting of the tank into user defined region
        (for variable meshing, for instance).

        If the user has defined no regions, this keeps the initial frame
        setup.
        """
        if not self.region_boundaries:
            return

        frame_segments = self.segments
        frame_segment_flags = self.segmentFlags

        self.segments = []
        self.segmentFlags = []

        new_segments = []
        new_segmentFlags = []
        self.region_column = {}
        corner_points = []

        for segment in frame_segments:
            segment_flag = frame_segment_flags[frame_segments.index(segment)]
            x_intersects = []
            segment_coords = self.vertices[segment[0]], \
                             self.vertices[segment[1]]
            for x_coord in self.region_boundaries:
                if min(segment_coords[0]) < x_coord < max(segment_coords[0]):
                    x_intersects += x_coord
                elif min(segment_coords[0]) == x_coord:
                    if max(segment_coords[0]) == x_coord:
                        continue
                    else:
                        corner_points += min(segment_coords)
                elif max(segment_coords[0]) == x_coord:
                    corner_points += max(segment_coords)
            start_point = self.vertices.index(min(segment_coords))
            end_point = self.vertices.index(max(segment_coords))
            for x in x_intersects:
                y = self._lineVertexInterpolation(x, segment)
                self.vertices += [[x, y], ]
                self.vertexFlags += [segment_flag, ]
                new_point = len(self.vertices) - 1
                new_segments += [[start_point, new_point], ]
                new_segmentFlags += [segment_flag, ]
                start_point = new_point
            new_segments += [[start_point, end_point], ]
            new_segmentFlags += [segment_flag, ]

        for x in self.region_boundaries:
            column_points = [vert for vert in self.vertices if vert[0] == x]
            column_points += corner_points
            column_points.sort(key=lambda v: self.vertices[v][1])
            y_n_points = column_points[1::2]
            y_p_points = column_points[2::2]
            if len(y_n_points) != len(y_p_points):
                raise IndexError('ERROR: Boundary at x coordinate '
                                 + str(x)
                                 + ' is confused. Shapes may be inappropriate')
            self.region_column = []
            for index in range(len(y_n_points)):
                y_n_index = self.vertices.index(y_n_points[index])
                y_p_index = self.vertices.index(y_p_points[index])
                new_segments += [[y_n_index, y_p_index], ]
                new_segmentFlags += [self.boundaryTags['empty'], ]
                self.region_column += [[y_n_index, y_p_index], ]

        self._addSegments(new_segments, new_segmentFlags, override = True)

    def _constructRegions(self):
        """Add in regions and holes after the final frame is made."""
        holes = []
        regions = []
        regionFlags = []

        if self.spongeLayers['x-']:
            regions += [[self.x0 - np.finfo(float).eps,
                         self.y1 - np.finfo(float).eps]]
            regionFlags += ['x-']
        if self.spongeLayers['x+']:
            regions += [[self.x1 + np.finfo(float).eps,
                         self.y1 - np.finfo(float).eps]]
            regionFlags += ['x+']
        if not self.region_boundaries:
            regions += [[self.x0 + np.finfo(float).eps,
                         self.y1 - np.finfo(float).eps]]
            # #[temp] to see if it fixes - a long term fix is needed later
            # regions[-1][1] = self.y1 - 0.1
            regionFlags += [1]
        else:
            region_num = 1
            for x in self.region_boundaries:
                for segment in self.region_column[x]:
                    halfway_height = [self.vertices[segment[0]][1]
                                      + self.vertices[segment[1]][1]
                                      ] * 0.5
                    x_n_side = [[x - np.finfo(float).eps,
                                 halfway_height]]
                    x_p_side = [[x + np.finfo(float).eps,
                                 halfway_height]]
                    regions += [x_n_side, x_p_side]
                    regionFlags += [region_num, region_num + 1]
                region_num += 1

        if self.floating_centers:
            for hole in self.floating_centers:
                holes += hole
            self.setHoles(holes)

        self.setRegions(regions, regionFlags)

    def constructShape(self):
        """Build the geometry (vertices,segments,regions) of the shape.
        """
        if not self.segments:
            self._constructInitialFrame()
        self._constructAdditionalFramework()
        self._constructRegions()

    def setDimensions(self, dim, new_coords=None):
        """Set dimension of the tank

        Parameters
        ----------
        dim: array_like
            Dimensions of the tank (without sponge layers)
        new_coords: Optional[array_like]
            If not None, will replace previous coords (or the calculated center
            of the rectangle, if from_0 == True) as the center.
        """
        self.dim = L, H = dim
        if self.from_0 is True:
            x, y = L / 2., H / 2.
        elif new_coords:
            x, y = new_coords
        else:
            x, y = self.coords
        self.coords = [x, y]
        self.x0, self.x1 = x - 0.5 * L, x + 0.5 * L
        self.y0, self.y1 = y - 0.5 * H, y + 0.5 * H

    def setHorizontalVariableMesh(self, region_boundaries):
        """Add region boundaries between potentially distinct refinements.

        Define points at which the mesh will be horizontally split, and
        may be given different regional constraints (or for any reason
        need to be different regions).

        Parameters
        ----------
        region_boundaries: array_like (floats)
            x-locations of additional regional boundaries inside
            the shape
        """
        for x in region_boundaries:
            if x < self.x0 or x > self.x1:
                raise ValueError('ERROR: boundary value' + str(x)
                                 + ' is not within the tank.')
        self.region_boundaries += region_boundaries
        self.region_boundaries.sort()
