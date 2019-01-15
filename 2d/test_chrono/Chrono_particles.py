import Chrono
import numpy as np


class Particles:
    def set_external_force(self,forces,torques):
        raise NotImplementedError

    def update_particle_states(self,dt):
        raise NotImplementedError

    def get_particle_states(self):
        raise NotImplementedError


class Chrono_particles(Particles):
    def __init__(self,
                timeStep,
                m_container_center,
                m_container_dims,
                m_gravity,
                nParticle,
                ball_center,
                ball_radius,
                ball_density):
        self.chmodel = Chrono.MBDModel(
                                        m_timeStep=timeStep,
                                        m_container_center=np.array(m_container_center,dtype="d"),
                                        m_container_dims=np.array(m_container_dims,dtype="d"),
                                        m_gravity=np.array(m_gravity,dtype="d"),
                                        nParticle = nParticle,
                                        ball_center=ball_center,
                                        ball_radius=ball_radius,
                                        ball_density=ball_density
                                        )
        self.force_from_fluid = np.zeros((nParticle,3),'d')
        self.torque_from_fluid= np.zeros((nParticle,3),'d')
        self.pos = np.zeros((nParticle,3),'d')
        self.vel= np.zeros((nParticle,3),'d')
        self.ang_vel= np.zeros((nParticle,3),'d')
        self.nParticle = nParticle

    def set_external_force(self,forces,torques):
        self.force_from_fluid = forces
        self.torque_from_fluid= torques

    def update_particle_states(self,dt):
        self.chmodel.step(self.force_from_fluid,self.torque_from_fluid,dt)
    
    def get_particle_states(self):
        self.chmodel.get_state(self.pos,self.vel,self.ang_vel)
        return self.pos,self.vel,self.ang_vel


def get_Chrono_particles():
    g_chrono = [0.0, -981, 0.0]

    L = [1.0, 1.0, 10.0]
    container_dim=[L[0],L[1],L[2]] #Dimensions of the container (Height/Width/depth)
    container_cent=[L[0]/2,L[1]/2,0.0] #Position of the center of the container"
    
    nParticle = 3
    ball_center = np.zeros((nParticle,3),'d')
    ball_radius = np.ones((nParticle,),'d')*0.1
    ball_density= np.ones((nParticle,),'d')*0.5

    ball_center[0,:] = 0.4,0.5,0.0
    ball_center[1,:] = 0.5,0.9,0.0
    ball_center[2,:] = 0.6,0.5,0.0

    dt_Chrono=0.001

    my_simulation = Chrono_particles(dt_Chrono,container_cent,container_dim,g_chrono,nParticle,
                                        ball_center,ball_radius,ball_density)
    return my_simulation


def test_particles():

    my_simulation = get_Chrono_particles()
    t = 0.0
    dt = 0.001
    forces = np.zeros((my_simulation.nParticle,3),'d')
    # forces[:,1] = 0.0157079632679*981##used to cancel gravity
    torques= np.zeros((my_simulation.nParticle,3),'d')
    with open('res.csv','w') as file:
        while t < 1.0:
            my_simulation.set_external_force(forces,torques)
            my_simulation.update_particle_states(dt)
            pos,vel,ang_vel = my_simulation.get_particle_states()
            print(pos,vel,ang_vel)
            t += dt

            ## Save data
            for i in range(my_simulation.nParticle):
                file.write('\t'.join([str(pos[i,j]) for j in range(3)]))
                file.write('\n')

if __name__ == "__main__":
    test_particles()