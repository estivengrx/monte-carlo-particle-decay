"""
Run analysis for muon decay Monte Carlo simulation.

This script generates and visualizes the electron energy spectrum from muon decay,
with and without radiative corrections, using the muonMonteCarlo class.
"""
from sys import path
from os.path import dirname, abspath

# Add project root to Python path
path.append(dirname(dirname(abspath(__file__))))
from src.python.muonMonteCarlo import muonMonteCarlo

def run_case(include_radiative, title, n_events):
    mc = muonMonteCarlo(n_events, include_radiative=include_radiative)
    mc.simulation(x_min=4.770585e-03, x_max=0.99999) # set x_min to be the lowest value that allows dGamma/dx to be non-negative, and x_max to be just below 1 to avoid the singularity at x=1
    mc.graphical_analysis(
        title=title,
        xtitle="Electron energy $E_e$ [MeV]",
        ytitle="Number of events"
    )
    print(mc.x_samples.min(), mc.x_samples.max())

def main():
    configs = [
        (True,  "Muon decay electron energy spectrum with radiative corrections", 100_000_000),
        (False, "Muon decay electron energy spectrum without radiative corrections", 100_000_000),
    ]
    for include_radiative, title, n_events in configs:
        run_case(include_radiative, title, n_events)

if __name__ == "__main__":
    main()