import numpy as np
from proteus import AuxiliaryVariables, Archiver, Comm, Profiling
cimport numpy
from proteus.mprans import BodyDynamics
from proteus import SpatialTools as st

cdef extern from "ChRigidBody.h":
    cdef cppclass ChVector:
        double x
        double y
        double z
    cdef cppclass ChQuaternion:
        double e0
        double e1
        double e2
        double e3
    cdef cppclass ChMatrix33:
        ChVector Get_A_Xaxis()
        ChVector Get_A_Yaxis()
        ChVector Get_A_Zaxis()

cdef extern from "ChRigidBody.h":
    cdef cppclass cppSystem:
        void DoStepDynamics(dt)
        void step(double dt)
    cppSystem * newSystem(double* gravity)

cdef extern from "ChRigidBody.h":
    cdef cppclass cppRigidBody:
        double mass
        ChVector pos
        ChQuaternion rot
        ChMatrix33 a
        # double* free_x
        # double* free_r
        void prestep(double* force, double* torque, double dt)
        void poststep()
        double hx(double* x, double dt)
        double hy(double* x, double dt)
        double hz(double* x, double dt)
    cppRigidBody * newRigidBody(cppSystem* system,
                                double* center,
                                double* rot,
                                double mass,
                                double* inertia,
                                double* free_x,
                                double* free_r)


cdef class RigidBody:
    cdef cppRigidBody * thisptr
    cdef public:
      object model
      object system
      object Shape
      numpy.ndarray barycneter0
      int nd, i_start, i_end
      double dt
      numpy.ndarray F
      numpy.ndarray M
      numpy.ndarray barycenter0
      # numpy.ndarray free_r
      # numpy.ndarray free_x
    def __cinit__(self,
                  shape,
                  System system,
                  numpy.ndarray center,
                  numpy.ndarray rot,
                  double mass,
                  numpy.ndarray inertia,
                  numpy.ndarray free_x,
                  numpy.ndarray free_r):
        self.system = system
        self.Shape = shape
        self.nd = shape.nd
        self.system.addBody(self)
        # rot = numpy.ndarray((4,))
        # rot[0] = np.sqrt(1+rot[0,0]+rot[1,1])/2.
        # rot[1] = 0#(rot[2,1]-rot[1,2])/(4*rot[0])
        # rot[2] = 0#(rot[0,2]-rot[2,0])/(4*rot[0])
        # rot[3] = (rot[1,0]-rot[0,1])/(4*rot[0])
        # rot[1] = (rot[2,1]-rot[1,2])/(4*rot[0])
        # rot[2] = (rot[0,2]-rot[2,0])/(4*rot[0])
        # rot[3] = (rot[1,0]-rot[0,1])/(4*rot[0])
        print(free_x)
        print(free_r)
        self.thisptr = newRigidBody(system.thisptr,
                                    <double*> center.data,
                                    <double*> rot.data,
                                    <double> mass,
                                    <double*> inertia.data,
                                    <double*> free_x.data,
                                    <double*> free_r.data)
        print('llala')
        if 'ChRigidBody' not in shape.auxiliaryVariables:
            shape.auxiliaryVariables['ChRigidBody'] = self

    def set_indices(self, i_start, i_end):
        self.i_start = i_start
        self.i_end = i_end

    def attachAuxiliaryVariables(self,avDict):
        pass


    def attachModel(self, model, ar):
        self.model = model
        return self

    def hx(self, numpy.ndarray x, double t):
        return self.thisptr.hx(<double*> x.data, t)

    def hy(self, numpy.ndarray x, double t):
        return self.thisptr.hy(<double*> x.data, t)

    def hz(self, numpy.ndarray x, double t):
        return self.thisptr.hz(<double*> x.data, t)

    # def setConstraintsDOF(self, numpy.ndarray free_x, numpy.ndarray free_r):
    #     """
    #     Sets constraints on the Shape (for moving bodies)

    #     Parameters
    #     ----------
    #     free_x: array_like
    #         Translational constraints.
    #     free_r: array_like
    #         Rotational constraints.
    #     """
    #     self.thisptr.free_x = <double*> free_x.data
    #     self.thisptr.free_r = <double*> free_r.data

    def setMass(self, mass):
        """
        Set mass of the shape.

        Parameters
        ----------
        mass: float
            mass of the body
        """
        self.thisptr.mass = <double> mass

    # def setInertiaTensor(self, It):
    #     """
    #     Set the inertia tensor of the shape

    #     Parameters
    #     ----------
    #     It: array_like, float
    #         Inertia tensor of the body (3x3 array in 3D, float in 2D)

    #     Notes
    #     -----
    #     The inertia tensor should not be already scaled with the mass of the
    #     shape.
    #     """
    #     self.thisptr.inertia = <double*> It

    def getPressureForces(self):
        """
        Gives the pressure forces applied on each segments/facets of the rigid
        body

        Returns
        -------
        F_p: array_like
            pressure forces (x, y, z) as provided by Proteus
        """
        i0, i1 = self.i_start, self.i_end
        F_p = self.model.levelModelList[-1].coefficients.netForces_p[i0:i1, :]
        print("FP : ", self.model.levelModelList[-1].coefficients.netForces_p)
        print("FP0:", F_p)
        F_t = np.sum(F_p, axis=0)
        return F_t

    def getShearForces(self):
        """
        Gives the shear forces applied on each segments/facets of the rigid
        body

        Returns
        -------
        F_v: array_like
            shear forces (x, y, z) as provided by Proteus
        """
        i0, i1 = self.i_start, self.i_end
        F_v = self.model.levelModelList[-1].coefficients.netForces_v[i0:i1, :]
        F_t = np.sum(F_v, axis=0)
        return F_t

    def getMoments(self):
        """
        Gives the moments applied on each segments/facets of the rigid body

        Returns
        -------
        M: array_like
            moments (x, y, z) as provided by Proteus
        """
        i0, i1 = self.i_start, self.i_end
        M = self.model.levelModelList[-1].coefficients.netMoments[i0:i1, :]
        M_t = np.sum(M, axis=0)
        # !!!!!!!!!!!! UPDATE BARYCENTER !!!!!!!!!!!!
        Fx, Fy, Fz = self.F
        rx, ry, rz = self.barycenter0-self.getPosition()
        Mp = np.array([ry*Fz-rz*Fy, -(rx*Fz-rz*Fx), (rx*Fy-ry*Fx)])
        M_t += Mp
        return M_t

    def getPosition(self):
        x = self.thisptr.pos.x
        y = self.thisptr.pos.y
        z = self.thisptr.pos.z
        return np.array([x, y, z])

    def getRotationQuaternion(self):
        e0 = self.thisptr.rot.e0
        e1 = self.thisptr.rot.e1
        e2 = self.thisptr.rot.e2
        e3 = self.thisptr.rot.e3
        return np.array([e0, e1, e2, e3])

    def getRotationMatrix(self):
        x0 = self.thisptr.a.Get_A_Xaxis().x
        x1 = self.thisptr.a.Get_A_Xaxis().y
        x2 = self.thisptr.a.Get_A_Xaxis().z
        y0 = self.thisptr.a.Get_A_Yaxis().x
        y1 = self.thisptr.a.Get_A_Yaxis().y
        y2 = self.thisptr.a.Get_A_Yaxis().z
        z0 = self.thisptr.a.Get_A_Zaxis().x
        z1 = self.thisptr.a.Get_A_Zaxis().y
        z2 = self.thisptr.a.Get_A_Zaxis().z
        matrix = np.array([x0, x1, x2],
                          [y0, y1, y2],
                          [z0, z1, z2])
        return matrix

    def step(self, numpy.ndarray F, numpy.ndarray M, double dt):
        self.thisptr.prestep(<double*> F.data, <double*> M.data, dt)

    def poststep(self):
        self.thisptr.poststep()

    def calculate_init(self):
        a = 1
        #

    def calculate(self):
        try:
            self.dt = self.model.levelModelList[-1].dt_last
        except:
            self.dt = self.system.dt_init
        self.F = self.getPressureForces()+self.getShearForces()
        self.M = self.getMoments()
        self.step(self.F, self.M, self.dt)


cdef class System:
    cdef cppSystem * thisptr
    cdef object model
    cdef object bodies
    cdef public double dt_init
    cdef double dt
    def __cinit__(self, numpy.ndarray gravity):
        self.thisptr = newSystem(<double*> gravity.data)
        self.bodies = []
        self.dt_init = 0.001

    def attachModel(self, model, ar):
        self.model = model
        return self
    def attachAuxiliaryVariables(self,avDict):
        pass
    def calculate_init(self):
        for body in self.bodies:
            body.barycenter0 = body.Shape.barycenter.copy()
        a = 1
        #

    def calculate(self):
        a = 1
        try:
            self.dt = self.model.levelModelList[-1].dt_last
        except:
            self.dt = self.dt_init
        dt = self.dt
        self.step(dt)
        for body in self.bodies:
            body.poststep()
        #

    def step(self, double dt):
        self.thisptr.step(<double> dt)

    def addBody(self, body):
        # !!!! TO CHANGE !!!!
        self.bodies += [body]
