# Code created using Copilot GPT-4.1 to quickly test the muonMonteCarlo class I created.

import unittest
import numpy as np
from sys import path
from os.path import dirname, abspath, join, exists

# Add project root to Python path
path.append(dirname(dirname(abspath(__file__))))

# Import using the full module path from the root
from src.python.muonMonteCarlo import muonMonteCarlo

class TestMuonMonteCarlo(unittest.TestCase):
    def setUp(self):
        """Set up a fresh instance before each test with fewer events for speed."""
        self.mc_sim = muonMonteCarlo(n_events=100)

    def test_initialization(self):
        """Test physical constants and instance variables are initialized properly."""
        self.assertEqual(self.mc_sim.n_events, 100)
        self.assertGreater(self.mc_sim.muon_mass, 0)
        self.assertGreater(self.mc_sim.electron_mass, 0)

    def test_differential_decay_rate(self):
        """Test if the internal physics equation returns expected shapes and valid floats."""
        x_vals = np.linspace(0.1, 0.9, 100)
        rate = self.mc_sim._dGamma_dx_differential_decay_rate(x_vals)
        self.assertEqual(rate.shape, x_vals.shape)
        # Ensure probabilities/rates are non-negative
        self.assertTrue(np.all(rate >= 0))
    
    def test_acceptance_rejection_method(self):
        """Test the acceptance-rejection method produces the correct number of samples and valid values."""
        x_vals = np.linspace(0.1, 0.9, 1000)
        samples = self.mc_sim.random_sampling_acceptance_rejection(x_vals, 50)
        self.assertEqual(len(samples), 50)
        self.assertTrue(np.all(samples >= 0))
        self.assertTrue(np.all(samples <= 1))

    def test_simulation_generates_samples(self):
        """Test the acceptance-rejection method outputs exactly n_events."""
        self.mc_sim.simulation(x_min=0.1, x_max=0.99)
        # Check arrays exist
        self.assertTrue(hasattr(self.mc_sim, 'x_samples'))
        self.assertTrue(hasattr(self.mc_sim, 'E_samples'))
        # Check correct sample sizes were drawn
        self.assertEqual(len(self.mc_sim.x_samples), 100)
        self.assertEqual(len(self.mc_sim.E_samples), 100)
        # Verify physical boundaries (Energy shouldn't exceed half the muon mass)
        max_E = self.mc_sim.muon_mass / 2
        self.assertTrue(np.all(self.mc_sim.E_samples <= max_E))
        self.assertTrue(np.all(self.mc_sim.x_samples <= 1.0))
    
    def test_image_output_is_results_plot(self):
        """Test if the plot is saved in the expected directory with the expected name."""
        from unittest.mock import patch
        with patch('matplotlib.pyplot.show'):
            self.mc_sim.simulation(x_min=0.1, x_max=0.99)
            self.mc_sim.graphical_analysis(
                title="electron_energy_spectrum_test",
                xtitle="Electron Energy $E_e$ [MeV]",
                ytitle="Probability Density"
            )
        expected_path = join(self.mc_sim.BASE_DIR, 'results', 'plots', 'electron_energy_spectrum_test.png')
        self.assertTrue(exists(expected_path), f"Expected plot not found at {expected_path}")

if __name__ == "__main__":
    # Run the tests
    unittest.main()