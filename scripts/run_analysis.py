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
    mc.simulation(x_min=0.001, x_max=0.999)
    mc.graphical_analysis(
        title=title,
        xtitle="Electron Energy $E_e$ [MeV]",
        ytitle="Probability Density"
    )

def main():
    configs = [
        (True,  "Muon decay electron energy spectrum with radiative corrections", 100_000),
        (False, "Muon decay electron energy spectrum without radiative corrections", 100_000),
    ]
    for include_radiative, title, n_events in configs:
        run_case(include_radiative, title, n_events)

if __name__ == "__main__":
    main()