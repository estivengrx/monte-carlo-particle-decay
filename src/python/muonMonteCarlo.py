"""
Reference paper for this code:
Czarnecki, A., Dowling, M., Garcia i Tormo, X., Marciano, W. J., & Szafron, R. (2014). Michel decay spectrum for a muon bound
to a nucleus. Physical Review D, 90 (9). https://doi.org/10.1103/physrevd.90.093002
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import scienceplots as sp

from tqdm import tqdm
from pathlib import Path
from particle import Particle
from scipy.constants import physical_constants
plt.style.use(['science', 'ieee'])
mpl.use('Qt5Agg')


class muonMonteCarlo:
    def __init__(self, n_events:int = 100_000, include_radiative: bool = True):
        """Add docstring and method outputs"""
        self.n_events = n_events
        self.include_radiative = include_radiative

        # Physical programming constants
        self.muon_mass = Particle.from_name("mu-").mass
        self.electron_mass = Particle.from_name("e-").mass
        self.Gf_fermi_constant = physical_constants['Fermi coupling constant'][0]
        self.alpha = physical_constants['fine-structure constant'][0]
        self.BASE_DIR = Path(__file__).resolve().parents[2] # Assuming this file is in src/python, we go up two levels to the project root

    # Internal Physics
    def _f_radiative_corrections(self, x: np.ndarray) -> np.ndarray:
        """
        Compute the radiative correction function f(x) for muon decay.

        This function implements the radiative corrections to the differential decay rate,
        as required for precision calculations in muon decay spectra.

        Parameters
        ----------
        x : np.ndarray
            Array of normalized electron energies (x = 2*E_e/m_mu), where 0 < x <= 1.

        Returns
        -------
        np.ndarray
            The radiative correction term evaluated at each value of x. The output array has the same shape as x.
        """
        from scipy.special import spence
        log_mass_ration_mu_e = np.log(self.muon_mass / self.electron_mass)
        
        term1 = ((5/(3*x**2)) + (16*x/3) + (4/x) + ((12 - 8*x) * np.log((1/x) - 1)) - 8) * log_mass_ration_mu_e
        
        term2 = (6 - 4*x) * (2 * spence(1 - x) - 2 * (np.log(x))**2 + np.log(x) + (np.log(1 - x) * (3 * np.log(x) - (1/x)- 1) - (np.pi**2)/3 - 2))
        
        term3 = ((1 - x) * (34 * x**2 + (5 - 34 * x**2 + 17 * x) * np.log(x) - 22*x)) / (3*x**2)
        
        term4 = 6 * (1 - x) * np.log(x)
        
        return term1 + term2 + term3 + term4
    
    def _dGamma_dx_differential_decay_rate(self, x: np.ndarray) -> np.ndarray:
        """
        Compute the differential decay rate dΓ/dx for muon decay as a function of the normalized electron energy x.

        Parameters
        ----------
        x : np.ndarray
            Array of normalized electron energies (x = 2*E_e/m_mu), where 0 < x <= 1.

        Returns
        -------
        np.ndarray
            The differential decay rate evaluated at each value of x. The output array has the same shape as x.
        """
        # Mask to ensure we only compute for valid x values (0 < x <= 1)
        x = np.asarray(x)
        result = np.zeros_like(x, dtype=float)
        mask = (x > 0) & (x <= 1)

        # Initial factor from the decay rate formula in the paper 
        initial_factor = (self.Gf_fermi_constant**2 * self.muon_mass**5) / (192 * (np.pi ** 3))
        
        # Radiative corrections term, which can be significant for precision calculations
        if self.include_radiative:
            radiative_term = (self.alpha / np.pi) * self._f_radiative_corrections(x[mask])
        else:
            radiative_term = 0

        # Results masked to valid x range
        result[mask] = initial_factor * (x[mask]**2) * (6 - 4*x[mask] + radiative_term)
        return result
    
    def random_sampling_acceptance_rejection(self, x_values: np.ndarray, n_samples: int) -> np.ndarray:
        """
        Perform random sampling using the acceptance-rejection method.
        Uses vectorized batches for performance and tqdm for progress tracking.

        Parameters
        ----------
        x_values : np.ndarray
            Array of x values to evaluate the differential decay rate on.
        n_samples : int
            The number of accepted samples to generate.

        Returns
        -------
        np.ndarray
            Accepted normalized electron energy samples.
        """
        max_rate = np.max(self._dGamma_dx_differential_decay_rate(x_values))
        batch_size = 500_000  # candidates generated per iteration
        accepted = []
        n_accepted = 0

        # tqdm tracks accepted samples, total = n_samples is the finish line
        with tqdm(total=n_samples,
                desc="  Sampling",
                unit=" events",
                colour="blue",
                dynamic_ncols=True) as pbar:

            while n_accepted < n_samples:
                x_rand = np.random.uniform(0, 1, batch_size)
                y_rand = np.random.uniform(0, max_rate, batch_size)
                mask = y_rand < self._dGamma_dx_differential_decay_rate(x_rand)

                batch_accepted = x_rand[mask]
                accepted.append(batch_accepted)

                new = len(batch_accepted)
                n_accepted += new
                pbar.update(new)   # advance bar by however many were accepted

        return np.concatenate(accepted)[:n_samples] # ensuring we return exactly n_samples, in case we accepted a few extra in the last batch 

    def simulation(self, x_min: float, x_max: float) -> None:
        """
        Run the Monte Carlo simulation for muon decay and generate electron energy samples.

        Parameters
        ----------
        x_min : float
            Minimum value of the normalized electron energy (x) to sample.
        x_max : float
            Maximum value of the normalized electron energy (x) to sample.

        Side Effects
        ------------
        Sets the following instance attributes:
        - self.x_samples : np.ndarray
            Array of sampled normalized electron energies.
        - self.E_samples : np.ndarray
            Array of sampled electron energies in MeV.
        """
        self.x_samples = self.random_sampling_acceptance_rejection(np.linspace(x_min, x_max, self.n_events), self.n_events)
        self.E_samples = self.x_samples * self.muon_mass / 2

        # Saving raw data for potential future analysis as npy files
        data_dir = self.BASE_DIR / 'data' / 'raw'
        data_dir.mkdir(parents=True, exist_ok=True)
        np.save(data_dir / f'muon_decay_samples_{"radiative" if self.include_radiative else "no_radiative"}_{self.n_events}_events.npy', self.E_samples)
    
    def error(self):
        """Add docstring and method outputs"""
        pass

    def comparation_real_data(self):
        """Add docstring and method outputs"""
        pass

    def graphical_analysis(self, title: str,
                                 xtitle: str,
                                 ytitle: str
                                 ) -> None:
        """
        Generate and display a scientific-quality histogram of the simulated electron energies,
        using pythons library SciencePlots.

        This method creates a histogram of the electron energy samples generated by the simulation,
        applies scientific plotting styles, and saves the resulting figure as a high-resolution PNG
        in the results/plots directory. The plot is also displayed interactively.

        Side Effects
        ------------
        - Saves the plot to '{BASE_DIR}/results/plots/electron_energy_spectrum.png'.
        - Displays the plot window.
        """
        fig, ax = plt.subplots(figsize=(6, 4))

        # Histogram with scientific style
        counts, bins, patches = ax.hist(
            self.E_samples,
            bins='fd', # Freedman-Diaconis rule for optimal binning
            density=True,
            linewidth=1,
            histtype='step',
            color='C0',
            label='Simulated spectrum'
        )

        # Axis labels
        ax.set_xlabel(xtitle)
        ax.set_ylabel(ytitle)

        # Title
        ax.set_title(title)

        # Minor ticks + grid
        ax.minorticks_on()
        ax.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.6)
        fig.tight_layout()
        ax.legend()

        # Figure saving
        plt.savefig(f'{self.BASE_DIR}/results/plots/{title.replace(" ", "_").lower()}.png', dpi=300)

    def export_data_to_geant4(self):
        """Add docstring and method outputs"""
        pass