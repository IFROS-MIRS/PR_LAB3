from math import pi
import numpy as np
import scipy
from MCLocalization import MCLocalization
from DifferentialDriveSimulatedRobot import *
from Pose3D import *
from ParticleFilter import pdf

class PF_3DOF(MCLocalization):
    def __init__(self, index, kSteps, robot, Map, particles, *args):
        super().__init__(index, kSteps, robot, particles, *args)

        self.dt = 0.1  # dt is the sampling time at which we iterate the DR
        self.wheelRadius = 0.1  # wheel radius
        self.wheelBase = 0.5  # wheel base
        self.robot.pulse_x_wheelTurns = 4096  # number of pulses per wheel turn

        self.M = Map # Save the map for the observation model (Map Based Localization)

    def GetInput(self):
        """
        Get the input for the motion model.

        :return: * **uk, Qk**. uk: input vector (:math:`u_k={}^B[\Delta L~\Delta R]^T`), Qk: covariance of the input noise

        **To be completed by the student**.
        """

        #
        #
        #return uk, Qk
    
    def GetMeasurements(self):
        """ 
        Get the measurements for the observation model.

        :return: * **zf, Rf**. zf: list of measurements, Rf: list of covariance matrices of the measurements

        **To be completed by the student**. 
        """

        return [], [] # Used to test the code without measurements!

    
    def SampleProcessModel(self, particle, u, Q):
        """ 
        Apply the process model to a single particle.

        **To be completed by the student**. 
        """

        #
        #
    
    
    def ObservationModel(self, particle, z, R):
        """ 
        Compute the measurement probability of a single particle with respect to a single measurement.
        
        **To be completed by the student**. 
        """

        #
        #
        #


if __name__ == '__main__':

    
    M = [np.array([[-40, 5]]).T,
           np.array([[-5, 40]]).T,
           np.array([[-5, 25]]).T,
           np.array([[-3, 50]]).T,
           np.array([[-20, 3]]).T,
           np.array([[40,-40]]).T]  # feature map. Position of 2 point features in the world frame.
    

    #Simulation:
    xs0 = np.zeros((6, 1))
    kSteps = 5000
    index = [IndexStruct("x", 0, None), IndexStruct("y", 1, None), IndexStruct("yaw", 2, 0)]
    robot = DifferentialDriveSimulatedRobot(xs0, M)  # instantiate the simulated robot object

    # Particle Filter
    x0 = Pose3D(np.zeros((3,1)))  # initial guess
    P0 = np.diag([1**2, 1**2, np.deg2rad(20)**2]) # Initial uncertainty, CHOSEN BY THE STUDENT
    n_particles = 100 # Number of particles, CHOSEN BY THE STUDENT

    # create array of n_particles particles distributed randomly around x0 with covariance P
    #Each particle musy be a Pose3D object!!!

    # TO BE COMPLETED BY THE STUDENT
    # particles = [
    # ...
    # for _ in range(n_particles)
    #]
    #

    particles = [
        Pose3D(x0 + np.random.normal(0, P0).diagonal().reshape((3, 1)))
        for _ in range(n_particles)
    ]

    # particles is a np.array of Pose3D objects
    usk=np.array([[0.5, 0.03]]).T
    pf = PF_3DOF(index, kSteps, robot, M, particles)
    pf.LocalizationLoop(x0, usk)

    exit(0)
