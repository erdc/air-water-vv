import os
import csv
import numpy as np
from proteus import AuxiliaryVariables, Archiver, Comm, Profiling
cimport numpy as cnp
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
    cdef cppclass ChBody:
        void SetRot(ChQuaternion &mrot)

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
        ChBody body
        void prestep(double* force, double* torque, double dt)
        void poststep()
        double hx(double* x, double dt)
        double hy(double* x, double dt)
        double hz(double* x, double dt)
        void addSpring(double stiffness,
                       double damping,
                       double* fairlead,
                       double* anchor,
                       double rest_length)
        void setRotation(double* quat)
        void setPosition(double* pos)
        void setConstraints(double* free_x, double* free_r)
        void setInertiaXX(double* inertia)

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
      cnp.ndarray barycneter0
      int nd, i_start, i_end
      double dt
      dict record_dict
      cnp.ndarray F
      cnp.ndarray M
      cnp.ndarray barycenter0
      cnp.ndarray rotation_init
      # cnp.ndarray free_r
      # cnp.ndarray free_x
    def __cinit__(self,
                  shape,
                  System system,
                  cnp.ndarray center,
                  cnp.ndarray rot,
                  double mass,
                  cnp.ndarray inertia,
                  cnp.ndarray free_x,
                  cnp.ndarray free_r):
        self.system = system
        self.Shape = shape
        self.nd = shape.nd
        self.system.addBody(self)
        self.record_dict = {}
        self.thisptr = newRigidBody(system.thisptr,
                                    <double*> center.data,
                                    <double*> rot.data,
                                    <double> mass,
                                    <double*> inertia.data,
                                    <double*> free_x.data,
                                    <double*> free_r.data)
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

    def hx(self, cnp.ndarray x, double t):
        return self.thisptr.hx(<double*> x.data, t)

    def hy(self, cnp.ndarray x, double t):
        return self.thisptr.hy(<double*> x.data, t)

    def hz(self, cnp.ndarray x, double t):
        return self.thisptr.hz(<double*> x.data, t)

    # def setConstraintsDOF(self, cnp.ndarray free_x, cnp.ndarray free_r):
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
    def addSpring(self, double stiffness, double damping, cnp.ndarray fairlead,
                  cnp.ndarray anchor, double rest_length):
        self.thisptr.addSpring(stiffness, damping, <double*> fairlead.data,
                               <double*> anchor.data, rest_length)

    def setPosition(self, cnp.ndarray position):
        self.thisptr.setPosition(<double*> position.data)
        
    def setRotation(self, cnp.ndarray quaternion):
        self.thisptr.setRotation(<double*> quaternion.data)

    def setConstraints(self, cnp.ndarray free_x, cnp.ndarray free_r):
        self.thisptr.setConstraints(<double*> free_x.data, <double*> free_r.data)

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

    def step(self, cnp.ndarray F, cnp.ndarray M, double dt):
        self.thisptr.prestep(<double*> F.data, <double*> M.data, dt)

    def poststep(self):
        self.thisptr.poststep()

    def calculate_init(self):
        # barycenter0 used for moment calculations
        self.barycenter0 = self.Shape.barycenter.copy()
        # self.thisptr.setRotation(<double*> self.rotation_init.data)
        #

    def calculate(self):
        try:
            self.dt = self.model.levelModelList[-1].dt_last
        except:
            self.dt = self.system.dt_init
        self.F = self.getPressureForces()+self.getShearForces()
        self.M = self.getMoments()
        self.step(self.F, self.M, self.dt)

    def getLastValues(self):
        self.acceleration
        self.acceleration_last
        self.velocity
        self.velocity_last
        self.rotation
        self.rotation_last
        self.F
        self.F_last
        self.M
        self.M_last
        self.ang_acceleration
        self.ang_acceleration_last
        self.ang_velocity
        self.ang_velocity_last
        self.inertia

    def setRecordValues(self, filename=None, all_values=False, pos=False,
                        rot=False, ang_disp=False, F=False, M=False,
                        inertia=False, vel=False, acc=False, ang_vel=False, ang_acc=False):
        """
        Sets the rigid body attributes that are to be recorded in a csv file
        during the simulation.
        Parameters
        ----------
        filename: Optional[string]
            Name of file, if not set, the file will be named as follows:
            'record_[shape.name].csv'
        all_values: bool
            Set to True to record all values listed below.
        time: bool
            Time of recorded row (default: True).
        pos: bool
            Position of body (default: False. Set to True to record).
        rot: bool
            Rotation of body (default: False. Set to True to record).
        ang_disp: array
            Angular displecement calculated during rigid body calculation step.
            Applied on the body in order to make it rotating.
        F: bool
            Forces applied on body (default: False. Set to True to record).
        M: bool
            Moments applied on body (default: False. Set to True to record).
        inertia: bool
            Inertia of body (default: False. Set to True to record).
        vel: bool
            Velocity of body (default: False. Set to True to record).
        acc: bool
            Acceleration of body (default: False. Set to True to record).
        ang_vel: array
            Angular velocity of body (default: False. Set to True to record).
        ang_acc: bool
            Angular acceleration of body (default: False. Set to True to record).
        Notes
        -----
        To add another value manually, add to dictionary self.record_dict:
        key: header of the column in .csv
        value: list of length 2: [variable name, index within variable]
                                                 (if no index, use None)
        e.g. self.record_dict['m']['mass', None]
        """
        if all_values is True:
            pos = rot = F = M = acc = vel = ang_acc = ang_vel = True
        if pos is True:
            self.record_dict['x'] = ['last_position', 0]
            self.record_dict['y'] = ['last_position', 1]
            self.record_dict['z'] = ['last_position', 2]
        if rot is True:
            self.record_dict['rx'] = ['last_rotation_euler', 0]
            self.record_dict['ry'] = ['last_rotation_euler', 1]
            self.record_dict['rz'] = ['last_rotation_euler', 2]
        if F is True:
            self.record_dict['Fx'] = ['F', 0]
            self.record_dict['Fy'] = ['F', 1]
            self.record_dict['Fz'] = ['F', 2]
            Fx = Fy = Fz = True
        if M is True:
            self.record_dict['Mx'] = ['M', 0]
            self.record_dict['My'] = ['M', 1]
            self.record_dict['Mz'] = ['M', 2]
        if acc is True:
            self.record_dict['ax'] = ['acceleration', 0]
            self.record_dict['ay'] = ['acceleration', 1]
            self.record_dict['az'] = ['acceleration', 2]
        if vel is True:
            self.record_dict['ux'] = ['velocity', 0]
            self.record_dict['uy'] = ['velocity', 1]
            self.record_dict['uz'] = ['velocity', 2]
        if ang_acc is True:
            self.record_dict['ang_ax'] = ['ang_acc', 0]
            self.record_dict['ang_ay'] = ['ang_acc', 1]
            self.record_dict['ang_az'] = ['ang_acc', 2]
        if ang_vel is True:
            self.record_dict['ang_ux'] = ['ang_vel', 0]
            self.record_dict['ang_uy'] = ['ang_vel', 1]
            self.record_dict['ang_uz'] = ['ang_vel', 2]
        if inertia is True:
            self.record_dict['intertia'] = ['inertia', None]

    def _recordValues(self):
        """
        Records values of rigid body attributes at each time step in a csv file.
        """
        comm = Comm.get()
        if comm.isMaster():
            record_file = os.path.join(Profiling.logDir, 'record_' + self.Shape.name + '.csv')
            t_last = self.model.stepController.t_model_last
            dt_last = self.model.levelModelList[-1].dt_last
            t = t_last-dt_last
            values_towrite = [t]
            if t == 0:
                headers = ['t']
                for key in self.record_dict:
                    headers += [key]
                with open(self.record_file, 'w') as csvfile:
                    writer = csv.writer(csvfile, delimiter=',')
                    writer.writerow(headers)
            for key, val in self.record_dict.iteritems():
                if val[1] is not None:
                    values_towrite += [getattr(self, val[0])[val[1]]]
                else:
                    values_towrite += [getattr(self, val[0])]
            with open(self.record_file, 'a') as csvfile:
                writer = csv.writer(csvfile, delimiter=',')
                writer.writerow(values_towrite)


cdef class System:
    cdef cppSystem * thisptr
    cdef object model
    cdef object bodies
    cdef public double dt_init
    cdef double dt
    def __cinit__(self, cnp.ndarray gravity):
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
            body.calculate_init()
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
    def calculate_init(self):
        for body in self.bodies:
            body.calculate_init()


    def step(self, double dt):
        self.thisptr.step(<double> dt)

    def addBody(self, body):
        # !!!! TO CHANGE !!!!
        self.bodies += [body]
