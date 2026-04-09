import numpy as np
import matplotlib.pyplot as plt
import scienceplots as sp

from particle import Particle
from scipy.constants import physical_constants
plt.style.use(['science', 'ieee'])

class muonMonteCarlo:
    def __init__(self, n_events:int = 100_000):
        """Add docstring and method outputs"""
        self.n_events = n_events

        # Physical programming constants
        self.muon_mass = Particle.from_name("mu-").mass
        self.electron_mass = Particle.from_name("e-").mass
        self.Gf_fermi_constant = physical_constants['Fermi coupling constant'][0]
        self.alpha = physical_constants['fine-structure constant'][0]

        # Mathematical constants
        self.pi = np.pi
    
    # Internal Physics
    def _f_radiative_corrections(self, x: np.ndarray) -> np.ndarray:
        """Add docstring and method outputs"""
        from scipy.special import spence
        log_mass_ration_mu_e = np.log(self.muon_mass / self.electron_mass)
        term1 = ((5/(3*x**2)) + (16*x/3) + (4/x) + ((12 - 8*x) * np.log((1/x) - 1)) - 8) * log_mass_ration_mu_e
        term2 = (6 - 4*x) * (2 * spence(1 - x) - 2 * (np.log(x))**2 + np.log(x) + (np.log(1 - x) * (3 * np.log(x) - (1/x)- 1) - (np.pi**2)/3 - 2))
        term3 = ((1 - x) * (34 * x**2 + (5 - 34 * x**2 + 17 * x) * np.log(x) - 22*x)) / (3*x**2)
        term4 = 6 * (1 - x) * np.log(x)
        return term1 + term2 + term3 + term4
    
    # Physics equations
    def _dGamma_dx_differential_decay_rate(self, x: np.ndarray, 
                                           include_radiative: bool = True):
        """Add docstring and method outputs"""

        initial_factor = (self.Gf_fermi_constant**2 * self.muon_mass**5) / (192 * (self.pi ** 3))
        if include_radiative:
            radiative_term = (self.alpha / self.pi) * self._f_radiative_corrections(x)
        else:
            radiative_term = 0

        return initial_factor * (x**2) * (6 - 4*x + radiative_term)

    def random_sampling(self):
        """Add docstring and method outputs"""
        pass

    def simulation(self):
        """Add docstring and method outputs"""
        pass

    def error(self):
        """Add docstring and method outputs"""
        pass

    def comparation_real_data(self):
        """Add docstring and method outputs"""
        pass

    def graphical_analysis(self):
        """Add docstring and method outputs"""
        pass

# Can this be generalized to create a monte carlo for any function?