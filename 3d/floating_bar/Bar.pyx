import numpy
cimport numpy
from proteus import AuxiliaryVariables, Archiver
from proteus.Profiling import  logEvent
import ode

cdef extern from "Bar.h":
    cdef cppclass cppChRigidBar:
        void step(double* force, double* torque, double dt)
        double hx(double* x, double dt)
        double hy(double* x, double dt)
        double hz(double* x, double dt)
    cppChRigidBar* newChRigidBar(double glass_thickness,
                                 double* bar_center,
                   		 double* bar_dim,
                                 double* L,
                                 double* gravity,
                                 double mass,
                                 double* inertia,
                                 double* free_x,
                                 double* free_r)

cdef class RigidBar:
    cdef cppChRigidBar* thisptr
    cdef object model
    def __cinit__(self,
                  double glass_thickness,
                  numpy.ndarray g=numpy.array((0,0,-9.8)),
                  numpy.ndarray bar_center=numpy.array((0.0,0.0,0.0)),
                  numpy.ndarray bar_dim=numpy.array((1.0,1.0,1.0)),
                  numpy.ndarray L=numpy.array((1.0,1.0,1.0)),
                  double mass=1.0,
                  numpy.ndarray inertia=numpy.array((1,1,1)),
                  numpy.ndarray free_x =numpy.array((1,1,1)),
                  numpy.ndarray free_r =numpy.array((1,1,1))):
        self.thisptr = newChRigidBar(glass_thickness,
                                     <double*> g.data,
                                     <double*> bar_center.data,
                                     <double*> bar_dim.data,
                                     <double*> L.data,
                                     mass,
                                     <double*> inertia.data,
                                     <double*> free_x.data,
                                     <double*> free_r.data)
    def attachModel(self,model,ar):
        self.model=model
        return self
    def get_u(self):
        return self.last_velocity[0]
    def get_v(self):
        return self.last_velocity[1]
    def get_w(self):
        return self.last_velocity[2]
    def calculate_init(self):
        self.calculate()
    def hx(self, numpy.ndarray x, double t):
        return self.thisptr.hx(<double*> x.data,t)
    def hy(self, numpy.ndarray x, double t):
        return self.thisptr.hy(<double*> x.data,t)
    def hz(self, numpy.ndarray x, double t):
        return self.thisptr.hz(<double*> x.data,t)
    def step(self,
             numpy.ndarray force,
             numpy.ndarray torque,
             double dt):
        self.thisptr.step(<double*> force.data,
                          <double*> torque.data,
                          dt)
    def calculate(self):
        import  numpy as np
        from numpy.linalg import inv
        import copy
        try:
            dt = self.model.levelModelList[-1].dt_last
        except:
            dt = 1.0e-8
        #cdef numpy.ndarray F;
        F = self.model.levelModelList[-1].coefficients.netForces_p[7,:] + self.model.levelModelList[-1].coefficients.netForces_v[7,:];
        #cdef numpy.ndarray M;
        M = self.model.levelModelList[-1].coefficients.netMoments[7,:]
        logEvent("x Force " +`self.model.stepController.t_model_last`+" "+`F[0]`)
        logEvent("y Force " +`self.model.stepController.t_model_last`+" "+`F[1]`)
        logEvent("z Force " +`self.model.stepController.t_model_last`+" "+`F[2]`)
        logEvent("x Moment " +`self.model.stepController.t_model_last`+" "+`M[0]`)
        logEvent("y Moment " +`self.model.stepController.t_model_last`+" "+`M[1]`)
        logEvent("z Moment " +`self.model.stepController.t_model_last`+" "+`M[2]`)
        logEvent("dt " +`dt`)
        scriptMotion=False
        linearOnly=False
        self.step( F,
                   M,
                   dt)
