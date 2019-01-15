import numpy
cimport numpy
from proteus import AuxiliaryVariables, Archiver
from proteus.Profiling import  logEvent
from cython.view cimport array as cvarray
# import ode
import numpy as np
# Import the C-level symbols of numpy
cimport numpy as np
from libcpp cimport bool

cdef extern from "Chrono.h":
    cdef cppclass cppMBDModel:
        void step(double* forces, double* torques, double dt)
        void set_state(double *_pos, double *_vel, double *_ang_vel)
        void get_state(double *_pos, double *_vel, double *_ang_vel)
        void save_state()
        void reset_state()
        void reset_timer()

        const int num_particles
        double *diam_
        double *pos_
        double *vel_
        double *angular_vel_
        void writeThisFrame()

    cppMBDModel* newMBDModel(double m_timeStep,
                             double* m_container_center,
                             double* m_container_dims,
                             double* m_gravity,
                             int nParticle,
                             double* ball_center,
                             double* ball_radius,
                             double* ball_density)


cdef class MBDModel:
    cdef cppMBDModel* thisptr
    cdef object model

    def __cinit__(self,
                  double m_timeStep,
                  np.ndarray m_container_center,
                  np.ndarray m_container_dims,
                  np.ndarray m_gravity,
                  int nParticle,
                  np.ndarray ball_center,
                  np.ndarray ball_radius,
                  np.ndarray ball_density,
                  ):
       
        self.thisptr =  newMBDModel(m_timeStep,
                                    <double*> m_container_center.data,
                                    <double*> m_container_dims.data,
                                    <double*> m_gravity.data,
                                    nParticle,
                                    <double*> ball_center.data,
                                    <double*> ball_radius.data,
                                    <double*> ball_density.data
                                     )


    def getNumParticles(self):
        return self.thisptr.num_particles

    def attachModel(self,model,ar):
        self.model=model
        return self

    def get_u(self):
        return 0
    def get_v(self):
        return 0
    def get_w(self):
        return 0
        
    def calculate_init(self):
        self.calculate()

    def getParticlesDiameter(self, i):
        return self.thisptr.diam_[i]

    def get_Angular_Vel(self, i):
        vel=np.array([ self.thisptr.angular_vel_[3*i+0],self.thisptr.angular_vel_[3*i+1], self.thisptr.angular_vel_[3*i+2] ])
        return vel

    def get_Vel(self, i):
        vel=np.array([ self.thisptr.vel_[3*i+0],self.thisptr.vel_[3*i+1], self.thisptr.vel_[3*i+2] ])
        return vel

    def get_Pos(self, i):
        pos=np.array([ self.thisptr.pos_[3*i+0],self.thisptr.pos_[3*i+1], self.thisptr.pos_[3*i+2] ])
        return pos

    def step(self,
             np.ndarray[double, ndim=2, mode="c"] force,
             np.ndarray[double, ndim=2, mode="c"] torque,
             double dt):
        self.thisptr.step(&force[0,0],
                          &torque[0,0],
                          dt)
    def set_state(self, 
                  np.ndarray[double, ndim=2, mode="c"] position,
                  np.ndarray[double, ndim=2, mode="c"] velocity,
                  np.ndarray[double, ndim=2, mode="c"] angular_velocity):
        self.thisptr.set_state(&position[0,0],
                               &velocity[0,0],
                               &angular_velocity[0,0])

    def get_state(self, 
                  np.ndarray[double, ndim=2, mode="c"] position,
                  np.ndarray[double, ndim=2, mode="c"] velocity,
                  np.ndarray[double, ndim=2, mode="c"] angular_velocity):
        self.thisptr.get_state(&position[0,0],
                               &velocity[0,0],
                               &angular_velocity[0,0])

    def calculate(self, 
                  np.ndarray[double, ndim=2, mode="c"] Forces,
                  np.ndarray[double, ndim=2, mode="c"] torques,
                  double dt):
        logEvent("Calling chrono with dt " +`dt`)
        self.step(Forces,torques,dt)

    def writeFrame(self):
        self.thisptr.writeThisFrame()
    def save_state(self):
        self.thisptr.save_state()
    def reset_state(self):
        self.thisptr.reset_state()
    def reset_timer(self):
        self.thisptr.reset_timer()
