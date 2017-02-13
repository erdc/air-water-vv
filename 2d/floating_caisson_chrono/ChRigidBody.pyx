# distutils: language = c++

import os
import csv
import numpy as np
from proteus import AuxiliaryVariables, Archiver, Comm, Profiling
cimport numpy as np
from proteus import SpatialTools as st
from libcpp.string cimport string
from libcpp cimport bool
from libcpp.memory cimport (shared_ptr,
                            make_shared)
from collections import OrderedDict


cdef extern from "ChRigidBody.h":
    cdef cppclass ChVector[double]:
        double x
        double y
        double z
    cdef cppclass ChQuaternion[double]:
        double e0
        double e1
        double e2
        double e3
    cdef cppclass ChMatrix33[double]:
        ChVector Get_A_Xaxis()
        ChVector Get_A_Yaxis()
        ChVector Get_A_Zaxis()
    cdef cppclass ChBody:
        void SetRot(ChQuaternion &rot)


cdef extern from "ChRigidBody.h":
    cdef cppclass cppSystem:
        void DoStepDynamics(dt)
        void step(double proteus_dt)
        void recordBodyList()
        void setChTimeStep(double dt)
        void setGravity(double* gravity)
        void setDirectory(string directory)
    cppSystem * newSystem(double* gravity)


cdef extern from "ChRigidBody.h":
    cdef cppclass cppRigidBody:
        double mass
        ChVector pos
        ChVector pos_last
        ChVector vel
        ChVector vel_last
        ChVector acc
        ChVector acc_last
        ChVector angvel
        ChVector angvel_last
        ChVector angacc
        ChVector angacc_last
        # ChVector inertia
        ChMatrix33 rotm
        ChMatrix33 rotm_last
        ChQuaternion rotq
        ChQuaternion rotq_last
        # double* free_x
        # double* free_r
        ChVector F
        ChVector F_last
        ChVector M
        ChVector M_last
        ChBody body
        cppRigidBody(cppSystem* system, double* position,
                     double* rotq, double mass, double* inertia,
                     double* free_x, double* free_r)
        void prestep(double* force, double* torque)
        void poststep()
        double hx(double* x, double dt)
        double hy(double* x, double dt)
        double hz(double* x, double dt)
        void addSpring(double stiffness,
                       double damping,
                       double* fairlead,
                       double* anchor,
                       double rest_length)
        void addPrismaticLinksWithSpring(double* pris1,
                                         double* pris2,
                                         double stiffness,
                                         double damping,
                                         double rest_length);
        void setRotation(double* quat)
        void setPosition(double* pos)
        void setConstraints(double* free_x, double* free_r)
        void setInertiaXX(double* inertia)
        void setName(string name)

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
      str record_file
      object model
      object system
      object Shape
      int nd, i_start, i_end
      double dt
      object record_dict
      object prescribed_motion_function
      np.ndarray position
      np.ndarray position_last
      np.ndarray F
      np.ndarray M
      np.ndarray F_last
      np.ndarray M_last
      np.ndarray acceleration
      np.ndarray acceleration_last
      np.ndarray velocity
      np.ndarray velocity_last
      np.ndarray ang_acceleration_last
      np.ndarray ang_acceleration
      np.ndarray ang_velocity_last
      np.ndarray ang_velocity
      np.ndarray barycenter0
      np.ndarray rotation_init
      np.ndarray rotm
      np.ndarray rotm_last
      np.ndarray rotq
      np.ndarray rotq_last
      np.ndarray F_prot
      np.ndarray M_prot
      np.ndarray F_prot_last
      np.ndarray M_prot_last
      # np.ndarray free_r
      # np.ndarray free_x
    def __cinit__(self,
                  shape,
                  System system,
                  np.ndarray center,
                  np.ndarray rot,
                  double mass,
                  np.ndarray inertia,
                  np.ndarray free_x,
                  np.ndarray free_r):
        self.system = system
        self.Shape = shape
        self.nd = shape.nd
        self.system.addBody(self)
        self.record_dict = OrderedDict()
        self.thisptr = newRigidBody(system.thisptr,
                                    <double*> center.data,
                                    <double*> rot.data,
                                    <double> mass,
                                    <double*> inertia.data,
                                    <double*> free_x.data,
                                    <double*> free_r.data)
        if 'ChRigidBody' not in shape.auxiliaryVariables:
            shape.auxiliaryVariables['ChRigidBody'] = self
        self.setName(shape.name)
        self.F_prot = np.zeros(3)
        self.M_prot = np.zeros(3)
        self.prescribed_motion_function = None

    def set_indices(self, i_start, i_end):
        self.i_start = i_start
        self.i_end = i_end

    def attachAuxiliaryVariables(self,avDict):
        pass

    def setInitialRot(self, rot):
        cdef np.ndarray zeros = np.zeros(3)
        self.rotation_init = rot
        self.thisptr.prestep(<double*> zeros.data,
                             <double*> zeros.data)
        if self.rotation_init is not None:
            Profiling.logEvent("$$$$$ SETTING ROT")
            self.thisptr.setRotation(<double*> self.rotation_init.data)
        self.thisptr.poststep()

    def hx(self, np.ndarray x, double t):
        return self.thisptr.hx(<double*> x.data, t)

    def hy(self, np.ndarray x, double t):
        return self.thisptr.hy(<double*> x.data, t)

    def hz(self, np.ndarray x, double t):
        return self.thisptr.hz(<double*> x.data, t)

    # def setConstraintsDOF(self, np.ndarray free_x, np.ndarray free_r):
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
    def addSpring(self, double stiffness, double damping, np.ndarray fairlead,
                  np.ndarray anchor, double rest_length):
        self.thisptr.addSpring(stiffness, damping, <double*> fairlead.data,
                               <double*> anchor.data, rest_length)

    def setPosition(self, np.ndarray position):
        self.thisptr.setPosition(<double*> position.data)

    def setRotation(self, np.ndarray quaternion):
        self.thisptr.setRotation(<double*> quaternion.data)

    def setConstraints(self, np.ndarray free_x, np.ndarray free_r):
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
        F_p = self.system.model.levelModelList[-1].coefficients.netForces_p[i0:i1, :]
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
        F_v = self.system.model.levelModelList[-1].coefficients.netForces_v[i0:i1, :]
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
        M = self.system.model.levelModelList[-1].coefficients.netMoments[i0:i1, :]
        M_t = np.sum(M, axis=0)
        # !!!!!!!!!!!! UPDATE BARYCENTER !!!!!!!!!!!!
        Fx, Fy, Fz = self.F_prot
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
        e0 = self.thisptr.rotq.e0
        e1 = self.thisptr.rotq.e1
        e2 = self.thisptr.rotq.e2
        e3 = self.thisptr.rotq.e3
        return np.array([e0, e1, e2, e3])

    def getRotationMatrix(self):
        x0 = self.thisptr.rotm.Get_A_Xaxis().x
        x1 = self.thisptr.rotm.Get_A_Xaxis().y
        x2 = self.thisptr.rotm.Get_A_Xaxis().z
        y0 = self.thisptr.rotm.Get_A_Yaxis().x
        y1 = self.thisptr.rotm.Get_A_Yaxis().y
        y2 = self.thisptr.rotm.Get_A_Yaxis().z
        z0 = self.thisptr.rotm.Get_A_Zaxis().x
        z1 = self.thisptr.rotm.Get_A_Zaxis().y
        z2 = self.thisptr.rotm.Get_A_Zaxis().z
        matrix = np.array([x0, x1, x2],
                          [y0, y1, y2],
                          [z0, z1, z2])
        return matrix

    def prestep(self):
        self.F_prot_last = np.array(self.F_prot)
        self.M_prot_last = np.array(self.M_prot)
        self.F_prot = self.getPressureForces()+self.getShearForces()
        self.M_prot = self.getMoments()
        try:
            dt = self.system.model.levelModelList[-1].dt_last
        except:
            dt = self.system.dt_init
        self.setExternalForces(self.F_prot, self.M_prot)

    def setExternalForces(self, np.ndarray forces, np.ndarray moments):
        self.thisptr.prestep(<double*> forces.data,
                             <double*> moments.data)

    def poststep(self):
        if self.prescribed_motion_function is not None:
            new_x = self.callPrescribedMotion(self.system.model.stepController.t_model_last)
            self.thisptr.setPosition(<double*> new_x.data)
        self.thisptr.poststep()
        self.getValues()
        comm = Comm.get()
        if comm.isMaster():
            self._recordValues()

    def calculate_init(self):
        # barycenter0 used for moment calculations
        self.barycenter0 = self.Shape.barycenter.copy()
        # get the initial values for F and M
        cdef np.ndarray zeros = np.zeros(3)
        self.setExternalForces(zeros, zeros)
        self.thisptr.poststep()
        self.getValues()
        # self.thisptr.setRotation(<double*> self.rotation_init.data)
        #

    def calculate(self):
        pass

    def setPrescribedMotion(self, function):
        self.prescribed_motion_function = function

    cdef np.ndarray callPrescribedMotion(self, double t):
        return self.prescribed_motion_function(t)

    def getValues(self):
        # position
        self.position = ChVector_to_npArray(self.thisptr.pos)
        self.position_last = ChVector_to_npArray(self.thisptr.pos_last)
        # rotation
        self.rotq = ChQuaternion_to_npArray(self.thisptr.rotq)
        self.rotq_last = ChQuaternion_to_npArray(self.thisptr.rotq_last)
        self.rotm = ChMatrix33_to_npArray(self.thisptr.rotm)
        self.rotm_last = ChMatrix33_to_npArray(self.thisptr.rotm_last)
        # acceleration
        self.acceleration = ChVector_to_npArray(self.thisptr.acc)
        self.acceleration_last = ChVector_to_npArray(self.thisptr.acc_last)
        # velocity
        self.velocity = ChVector_to_npArray(self.thisptr.vel)
        self.velocity_last = ChVector_to_npArray(self.thisptr.vel_last)
        #angular acceleration
        self.ang_acceleration = ChVector_to_npArray(self.thisptr.angacc)
        self.ang_acceleration_last = ChVector_to_npArray(self.thisptr.angacc_last)
        # angular velocity
        self.ang_velocity = ChVector_to_npArray(self.thisptr.angvel)
        self.ang_velocity_last = ChVector_to_npArray(self.thisptr.angvel_last)
        # force
        self.F = ChVector_to_npArray(self.thisptr.F)
        self.F_last = ChVector_to_npArray(self.thisptr.F_last)
        # moment
        self.M = ChVector_to_npArray(self.thisptr.M)
        self.M_last = ChVector_to_npArray(self.thisptr.M_last)
        # self.M_last
        # # self.inertia = ChVector_to_npArray(self.thisptr.)



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
            self.record_dict['x'] = ['position_last', 0]
            self.record_dict['y'] = ['position_last', 1]
            self.record_dict['z'] = ['position_last', 2]
        # if rot is True:
        #     self.record_dict['rx'] = ['last_rotation_euler', 0]
        #     self.record_dict['ry'] = ['last_rotation_euler', 1]
        #     self.record_dict['rz'] = ['last_rotation_euler', 2]
        if rot is True:
            self.record_dict['rotq_e0'] = ['rotq_last', 0]
            self.record_dict['rotq_e1'] = ['rotq_last', 1]
            self.record_dict['rotq_e2'] = ['rotq_last', 2]
            self.record_dict['rotq_e3'] = ['rotq_last', 3]
            # self.record_dict['rotm_a11'] = ['rotm_last', (0,0)]
            # self.record_dict['rotm_a12'] = ['rotm_last', (0,1)]
            # self.record_dict['rotm_a13'] = ['rotm_last', (0,2)]
            # self.record_dict['rotm_a21'] = ['rotm_last', (1,0)]
            # self.record_dict['rotm_a22'] = ['rotm_last', (1,1)]
            # self.record_dict['rotm_a23'] = ['rotm_last', (1,2)]
            # self.record_dict['rotm_a31'] = ['rotm_last', (2,0)]
            # self.record_dict['rotm_a32'] = ['rotm_last', (2,1)]
            # self.record_dict['rotm_a33'] = ['rotm_last', (2,2)]
        if F is True:
            self.record_dict['Fx'] = ['F', 0]
            self.record_dict['Fy'] = ['F', 1]
            self.record_dict['Fz'] = ['F', 2]
            self.record_dict['Fx_prot'] = ['F_prot', 0]
            self.record_dict['Fy_prot'] = ['F_prot', 1]
            self.record_dict['Fz_prot'] = ['F_prot', 2]
            Fx = Fy = Fz = True
        if M is True:
            self.record_dict['Mx'] = ['M', 0]
            self.record_dict['My'] = ['M', 1]
            self.record_dict['Mz'] = ['M', 2]
            self.record_dict['Mx_prot'] = ['M_prot', 0]
            self.record_dict['My_prot'] = ['M_prot', 1]
            self.record_dict['Mz_prot'] = ['M_prot', 2]
        if acc is True:
            self.record_dict['ax'] = ['acceleration_last', 0]
            self.record_dict['ay'] = ['acceleration_last', 1]
            self.record_dict['az'] = ['acceleration_last', 2]
        if vel is True:
            self.record_dict['ux'] = ['velocity_last', 0]
            self.record_dict['uy'] = ['velocity_last', 1]
            self.record_dict['uz'] = ['velocity_last', 2]
        if ang_acc is True:
            self.record_dict['ang_ax'] = ['ang_acceleration_last', 0]
            self.record_dict['ang_ay'] = ['ang_acceleration_last', 1]
            self.record_dict['ang_az'] = ['ang_acceleration_last', 2]
        if ang_vel is True:
            self.record_dict['ang_ux'] = ['ang_velocity_last', 0]
            self.record_dict['ang_uy'] = ['ang_velocity_last', 1]
            self.record_dict['ang_uz'] = ['ang_velocity_last', 2]
        if inertia is True:
            self.record_dict['intertia'] = ['inertia', None]

    def _recordValues(self):
        """
        Records values of rigid body attributes at each time step in a csv file.
        """
        self.record_file = os.path.join(Profiling.logDir, 'record_' + self.Shape.name + '.csv')
        t_last = self.system.model.stepController.t_model_last
        dt_last = self.system.model.levelModelList[-1].dt_last
        t = t_last-dt_last
        t_prot = Profiling.time()-Profiling.startTime
        values_towrite = [t, t_prot]
        if t == 0:
            headers = ['t', 't_prot']
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

    def addPrismaticLinksWithSpring(self, np.ndarray pris1,
                                    np.ndarray pris2, double stiffness, double damping,
                                    double rest_length):
        """
        fairlead: barycenter coords
        pris: absolute coords

        pris1-------fairlead(barycenter)
        |
        |
        |
        |
        pris2
        """
        self.thisptr.addPrismaticLinksWithSpring(<double*> pris1.data, 
                                                 <double*> pris2.data,
                                                 stiffness,
                                                 damping,
                                                 rest_length)

    def setName(self, string name):
        self.thisptr.setName(name)


cdef class System:
    cdef cppSystem * thisptr
    cdef public object model
    cdef object bodies
    cdef public double dt_init
    cdef double proteus_dt
    cdef double chrono_dt
    cdef string directory
    def __cinit__(self, np.ndarray gravity):
        self.thisptr = newSystem(<double*> gravity.data)
        self.bodies = []
        self.dt_init = 0.001

    def attachModel(self, model, ar):
        self.model = model
        return self

    def attachAuxiliaryVariables(self,avDict):
        pass

    def calculate(self):
        try:
            self.proteus_dt = self.model.levelModelList[-1].dt_last
        except:
            self.proteus_dt = self.dt_init
        for body in self.bodies:
            body.prestep()
        self.step(self.proteus_dt)
        for body in self.bodies:
            body.poststep()
        #self.recordBodyList()

    def calculate_init(self):
        self.directory = str(Profiling.logDir)+'/'
        self.thisptr.setDirectory(self.directory)
        for body in self.bodies:
            body.calculate_init()
        #self.recordBodyList()

    def setTimeStep(self, double dt):
        """Sets time step for Chrono solver.
        Calculations in Chrono will use this time step within the
        Proteus time step (if bigger)

        Parameters
        ----------
        dt: float
            time step
        """
        self.chrono_dt = dt
        self.thisptr.setChTimeStep(dt)

    def setGravity(self, np.ndarray gravity):
        self.thisptr.setGravity(<double*> gravity.data)

    def step(self, double dt):
        self.thisptr.step(<double> dt)

    def addBody(self, body):
        # !!!! TO CHANGE !!!!
        self.bodies += [body]

    def recordBodyList(self):
        comm = Comm.get()
        if comm.isMaster():
            self.thisptr.recordBodyList()
# ctypedef np.ndarray vecarray(ChVector)

# ctypedef np.ndarray (*ChVector_to_npArray) (ChVector)
cdef np.ndarray ChVector_to_npArray(ChVector &myvec):
    return np.array([myvec.x, myvec.y, myvec.z])

cdef np.ndarray ChQuaternion_to_npArray(ChQuaternion &quat):
    return np.array([quat.e0, quat.e1, quat.e2, quat.e3])

cdef np.ndarray ChMatrix33_to_npArray(ChMatrix33 &mat):
    return np.array([[mat.Get_A_Xaxis().x, mat.Get_A_Xaxis().y, mat.Get_A_Xaxis().z],
                     [mat.Get_A_Yaxis().x, mat.Get_A_Yaxis().y, mat.Get_A_Yaxis().z],
                     [mat.Get_A_Zaxis().x, mat.Get_A_Zaxis().y, mat.Get_A_Zaxis().z]])

#def testx():
#    cdef ChSystem system = ChSystem()
#    cdef ChBody bod = ChBody()
#    cdef ChVector oo = ChVector[double](2.,3.,4.)
#    bod.SetPos_dt(oo)
#    cdef ChVector& gg = bod.GetPos_dt()
#    print(gg.x, gg.y, gg.z)
