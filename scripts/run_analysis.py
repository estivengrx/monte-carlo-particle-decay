"""
Run analysis for muon decay Monte Carlo simulation.

This script generates and visualizes the electron energy spectrum from muon decay,
with and without radiative corrections, using the muonMonteCarlo class.
"""
from sys import path
from os.path import dirname, abspath
from numpy import concatenate

# Add project root to Python path
path.append(dirname(dirname(abspath(__file__))))
from src.python.muonMonteCarlo import muonMonteCarlo

def run_case(include_radiative, title, n_events, n_runs=100):
    all_samples = []

    for _ in range(n_runs):
        mc = muonMonteCarlo(n_events, include_radiative=include_radiative)
        mc.simulation(x_min=4.770585e-03, x_max=0.99999)
        all_samples.append(mc.E_samples)

    # Concatenate all runs → final MC sample
    all_samples = concatenate(all_samples)

    # Plot final histogram
    mc.E_samples = all_samples  # reuse plotting method
    mc.graphical_analysis(
        title=title,
        xtitle="Electron energy $E_e$ [MeV]",
        ytitle="Number of events"
    )

    print("Final MC range:", all_samples.min(), all_samples.max())


def main():
    configs = [
        (True,  "Muon decay spectrum with radiative corrections", 100_000, 100),
        (False, "Muon decay spectrum without radiative corrections", 100_000, 100),
    ]

    for include_radiative, title, n_events, n_runs in configs:
        run_case(include_radiative, title, n_events, n_runs)


if __name__ == "__main__":
    main()
