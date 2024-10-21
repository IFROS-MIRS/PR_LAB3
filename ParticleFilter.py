import numpy as np
import matplotlib.pyplot as plt
import math
import random

def pdf(mean, sigma, x):
    """Compute the PDF for a normal distribution."""
    return 1 / (sigma * math.sqrt(2 * math.pi)) * math.exp(- (x-mean)**2 / (2 * sigma**2))


class ParticleFilter:
    """
    Particle Filter Localization.

    This class implements basic plotting and logging functionality for the Particle Filter,
    as well as the interface for the child classes to implement.

    A particle filter is a Monte Carlo algorithm that approximates the posterior distribution of the robot
    by a set of weighted particles.  Note that the "weight" (which is a terrible term) is simply the 
    probability of the particle being correct. Therefore, each particle is an estimate, and each estimate 
    has some probability of being correct.

    """
    def __init__(self, particles, *args):
        """
        Constructor of the Particle Filter class.

        :param index: Logging index structure (:class:`Index`)
        :param kSteps: Number of time steps to simulate
        :param robot: Simulation robot object (:class:`Robot`)
        :param particles: initial particles as a list of Pose objects (or at least a list of numpy arrays)
        :param args: Rest of arguments to be passed to the parent constructor
        """
        self.particles = particles
        self.particle_weights = np.ones(len(particles)) / len(particles) # evenly distributed weights
        self.resampling_function = self.StochasticUniversalResampling # choose resampling method
        self.resampling_threshold = 10 # choose a resampling threshold


    def SampleProcessModel(self, p, u, Q):
        """
        Process model of a single particle given an input
        
        :param p: a single particle at time k
        :param u: input variable at time k
        :param Q: Covariance matrix associated to the input at time k

        :return the particle state at time k+1, after apllying the model
        """
        
        print("SampleProcessModel must be implemented in a child class")
        pass
    
    def ObservationModel(self, p, z, R):
        """
        Compute the measurement probability of a single particle with respect to a single measurement.

        :param p: a single particle
        :param z: a single measurement
        :param R: the covariance matrix associated to z

        :return: the measurement probability (or weight) of the particle p measuring z
        """
        print("ObservationModel must be implemented in a child class")
        pass


    def Prediction(self, u, Q):
        """
        Predict the next state of the system based on a given motion model.

        This function updates the state of each particle by predicting its next state using a motion model.

        :param u: input vector
        :param Q: the covariance matrix associated with the input vector
        :return: None

        **To be completed by the student**
        """

        # For each particle
        # --- Apply the process model to the particle



    def Update(self, z, R):    
        """
        Update the particle weights based on sensor measurements and perform resampling.

        This function adjusts the weights of particles based on how well they match the sensor measurements.
       
        The updated weights reflect the likelihood of each particle being the true state of the system given
        the sensor measurements.

        The resulting weights must be normalized

        After updating the weights, the function may perform resampling to ensure that particles with higher
        weights are more likely to be selected, maintaining diversity and preventing particle degeneracy.
        
        :param z: measurement vector
        :param R: the covariance matrix associated with the measurement vector

        :return: None

        **To be completed by the student**
        """

        # Update Weights of each particle based on the measurements
        #
        #
        #
        #
        #

        # Normalize weights
        #
        #

        # Resample if necessary
        # n_eff = ...
        #
        #


    def Resample(self) -> None:
        """
        Resample the particles based on their weights to ensure diversity and prevent particle degeneracy.

        This function implements the resampling step of a particle filter algorithm. It uses the weights
        assigned to each particle to determine their likelihood of being selected. Particles with higher weights
        are more likely to be selected, while those with lower weights have a lower chance.

        The resampling process helps to maintain a diverse set of particles that better represents the underlying
        probability distribution of the system state. 

        After resampling, the attributes 'particles' and 'weights' of the ParticleFilter instance are updated
        to reflect the new set of particles and their corresponding weights.

        This method calls the resampling function set as default, mainly either RouletteWheelResampling or StochasticUniversalResampling
       
        :return: None
        """
        return self.resampling_function()
    
    # Resampling methods

    def RouletteWheelResampling(self):
        ''' This method is the Roulette Wheel version of resampling'''

        '''**To be completed by the student**'''

        #
        #
        #...
        #
        #
        # self.particles = new_particles
        # self.particle_weights = ...

    def StochasticUniversalResampling(self):
        ''' This method is the Stochastic Universal version of resampling'''

        '''**To be completed by the student**'''

        #
        #
        #...
        #
        #
        # self.particles = new_particles
        # self.particle_weights = ...

    # Helper functions for plotting and logging

    def get_mean_particle(self):
        """
        Calculate the mean particle based on the current set of particles and their weights.
        :return: mean particle
        """
        # Weighted mean
        return np.average(self.particles, axis=0, weights=self.particle_weights)
    
    def get_best_particle(self):
        """
        Calculate the best particle based on the current set of particles and their weights.
        :return: best particle
        """
        # Maximum weight
        return self.particles[np.argmax(self.particle_weights)]

