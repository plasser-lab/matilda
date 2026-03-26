"""
version 1.0.0
author: Giulia Woelfle-Conway
usage: Routines for Franck-Condon simulations
"""

import math
from dataclasses import dataclass
from typing import Optional

import numpy
import matplotlib.pyplot as plt

from matilda import units, error_handler

class HRData:
    def __init__(self, fc_options, hr_data_file):
        self.hr_data_file = hr_data_file
        self.fc_options = fc_options

        self.freqs_au = []
        self.S_factors = []

        self.sig_modes = []
        self.FC_modes = []

        self.stick_energies = []
        self.stick_intensities = []

        self.sigma = None

    def read_hr_data(self):
        """
        Reads the data from hr.txt
        """
        for line in open(self.hr_data_file, 'r').readlines()[1:]:  # Skip headers
            parts = line.strip().split()
            if len(parts) < 3:              #ensures the file is correct
                continue
            freq_cm = float(parts[1])
            S = float(parts[2])

            freq_au = freq_cm / units.energy['rcm']
            self.freqs_au.append(freq_au)           #data for the containers previously made
            self.S_factors.append(S)

        if not self.freqs_au or not self.S_factors:
            raise error_handler.MsgError("Error: No valid data found in Huang-Rhys data file")

    def classify_modes(self):
        """
        Separates into Franck-Condon and broadening (sigma) modes
        """
        w_min_au = self.fc_options["w_min"] / units.energy['rcm']

        for imode in range(len(self.freqs_au)):
            omega = self.freqs_au[imode]
            S     = self.S_factors[imode]
            if S >= self.fc_options["S_min"] and omega >= w_min_au:
                self.FC_modes += [(S, omega)]            #runs through the list of values and assigns them according to the criteria
            else:
                self.sig_modes += [(S, omega)]

        if len(self.FC_modes) == 0:
            raise error_handler.MsgError("No modes above the given S_min and w_min for FC progression")

        self.sigma = math.sqrt(sum(0.5 * (omega**2) * S for S, omega in self.sig_modes))

    def prep(self):
        self.read_hr_data()
        self.classify_modes()
        #self.E_0_au = self.fc_options["dEH"]
    
        print(f'Using modes with omega > {self.fc_options["w_min"]:.4f} cm^-1 and S > {self.fc_options["S_min"]:.4f} as FC-active modes.')
        print(f"Number of active modes : {len(self.FC_modes)}")
        print(f"Number of sigma modes  : {len(self.sig_modes)}")
        print(f"Computed gaussian sigma: {self.sigma * units.energy['rcm']:.4f} cm^-1")
        print(f"Adiabatic energy difference: {self.fc_options['dEH'] * units.energy['rcm']:.4f} cm^-1")

    def compute_k_max(self):
        self.k_max = int(max(S + 5.0 * math.sqrt(S) for S, _ in self.FC_modes))

    def make_energy_grid(self):
        """
        Preparation of entire energy space
        """    
        fc_width = sum(S * omega for S, omega in self.FC_modes)

        if self.fc_options["abs_emi"] == 2:   
            #E_min = self.E_0_cm - 10000.0
            #E_max = self.E_0_cm + 10000.0
            E_min = self.fc_options["dEH"] - fc_width - 50.0 * self.sigma
            E_max = self.fc_options["dEH"] + 10.0 * self.sigma
        else:
            E_min = self.fc_options["dEH"] - 5.0 * self.sigma
            E_max = self.fc_options["dEH"] + fc_width + 50.0 * self.sigma

        self.energy_grid = numpy.linspace(E_min, E_max, self.fc_options['npoints'])

    def gaussian(self, x, x0):
        """
        Gaussian line shape in cm^-1
        """
        return numpy.exp(-0.5 * ((x - x0)/self.sigma)**2) / (self.sigma * numpy.sqrt(2 * numpy.pi))
    
    def build_spectrum(self):
        """
        Recursively enumerates vibronic quanta up to k_max per mode and places Gaussians.
        """
        self.spectrum = numpy.zeros_like(self.energy_grid)                        #accumulator for spectrum values computed later, one value per point in energy_grid
        S_total = sum(S for S, _ in self.FC_modes)
        prefactor = math.exp(-S_total)

        def loop_over_modes(idx, E_shift, intensity):                   #idx is the mode index currently being processed
            if idx == len(self.FC_modes):                                       # All modes handled -> place peak
                
                E_transition = self.fc_options["dEH"] + E_shift
                self.spectrum[:] += intensity * self.gaussian(self.energy_grid, E_transition)       #adds a gaussian peak centred at e_transition on the spectrum array and multiplying by intensity scales that peak

                self.stick_energies.append(E_transition)                     #records the discrete line at e_transition with weight intensity
                self.stick_intensities.append(intensity)
                return                                                  #stops when all modes are handled
                
            S, omega = self.FC_modes[idx]                                       #for current mode unpacks s and omega
            for k in range(self.k_max + 1):                                  #loops over no of quanta k in this mode from 0 up to k_max inclusive
                new_intensity = intensity * (S**k) / math.factorial(k)  #poisson-like weight factor is multiplied by intensity accumulated so far from previous modes. updates intenisty for choosing k quanta in this mode
                if self.fc_options["abs_emi"] == 2:                                  #updates energy shift depending on emi or abs
                    new_Eshift = E_shift - k * omega
                else:
                    new_Eshift = E_shift + k * omega
                loop_over_modes(idx + 1, new_Eshift, new_intensity)     #processes the next mode with updated e_shift and intensity. over full recursion, this enumerates all comb.s of k values for all modes

        loop_over_modes(0, 0.0, prefactor)                              #starts recursion at first mode, no energy shift yet and initial intensity at the prefactor

    def plot_postrun(self, plot_sticks: bool = True):
        """
        Prints out the results and saves data files
        """
        def write_table(filename: str, header: str, cols):
            with open(filename, "w") as fh:
                fh.write(header + "\n")
                for row in zip(*cols):
                    fh.write(" ".join(f"{float(v):15.8f}" for v in row) + "\n")
        
        def plot_and_save(x, y, xlabel, ylabel, title, filename, color, sticks=None, xlim=None,):
            plt.figure(figsize=(8, 5))
            plt.plot(x, y, color=color, label="Broadened spectrum")
            if sticks is not None:
                sx, sy = sticks
                plt.vlines(sx, 0, sy, color="blue", linewidth=1, label="Stick spectrum")
                plt.legend()
            if xlim is not None:
                plt.xlim(*xlim)
            plt.xlabel(xlabel)
            plt.ylabel(ylabel)
            plt.title(title)
            plt.grid(True, linestyle="--", alpha=0.5)
            plt.tight_layout()
            plt.savefig(filename, dpi=300)
            print(f"Plot saved as '{filename}'")

        wavenumber_grid = self.energy_grid * units.energy['rcm']
        wavelength_grid = units.energy['nm'] / self.energy_grid
        electronvolts_grid = self.energy_grid * units.energy['eV']

        stick_intensities_normalised = numpy.array(self.stick_intensities) / max(self.stick_intensities)

        stick_wavenumbers = numpy.array(self.stick_energies) * units.energy['rcm']
        stick_wavelengths = numpy.array(self.stick_energies) * units.energy['nm']
        stick_evolts = numpy.array(self.stick_energies) * units.energy['eV']

        xlim_cm = (float(wavenumber_grid.min()), float(wavenumber_grid.max()))
        xlim_nm = (float(wavelength_grid.min()), float(wavelength_grid.max()))
        xlim_ev = (float(electronvolts_grid.min()), float(electronvolts_grid.max()))

        #Emission
        if self.fc_options["abs_emi"] == 2:
            A_si = 10000 * 2 * math.pi * units.SIconstants['e']**2 / (units.SIconstants['me'] * units.SIconstants['eps0'] * units.SIconstants['c']) * self.fc_options["f"] * (self.fc_options["dEH"] * units.energy['rcm'])**2
            Lamda = 1e9 / A_si

            ein_coeff_emi = 2 * self.fc_options["dEH"]**2 * self.fc_options["f"] / units.constants['c0']**3 / units.time['s']
            lifetime = 1e9 / ein_coeff_emi

            self.spectrum /= self.spectrum.max()

            print("\nLineshapes are ignored.")
            #print(f"Einstein coefficient of spontaneous emission (A): {A_si:.6e} s^-1")

            print(f"Excited state lifetime: {Lamda:.6e} ns")
            print(f"Excited state lifetime au: {lifetime:.6e} ns")

            write_table("vibronic_emission_data.txt", "Wavenumber(cm^-1)   Wavelength(nm)  Electronvolts(eV)   Intensity(normalised)", [wavenumber_grid, wavelength_grid, electronvolts_grid, self.spectrum],)
            print("\nEmission spectrum data saved to 'vibronic_emission_data.txt'")

            if plot_sticks:
                write_table("vibronic_emi_stick_data.txt", "Wavenumber(cm^-1)  Wavelength(nm)   Electronvolts(eV)   Intensity(normalised)", [stick_wavenumbers, stick_wavelengths, stick_evolts, self.stick_intensities],)
                print("Stick spectrum data saved to vibronic_emi_stick_data.txt'")
            
            sticks_cm = (stick_wavenumbers, stick_intensities_normalised) if plot_sticks else None
            sticks_nm = (stick_wavelengths, stick_intensities_normalised) if plot_sticks else None
            sticks_ev = (stick_evolts, stick_intensities_normalised) if plot_sticks else None

            plot_and_save(wavenumber_grid, self.spectrum, "Wavenumber (cm$^{-1}$)", "Normalised Intensity", "Simulated Vibronic Emission Spectrum", "vibronic_emission.png", color="darkred", sticks=sticks_cm, xlim=xlim_cm,)
            plot_and_save(wavelength_grid, self.spectrum, "Wavelength (nm)", "Normalised Intensity", "Simulated Vibronic Emission Spectrum (Wavelength)", "vibronic_emission_nm.png", color="orange", sticks=sticks_nm, xlim=xlim_nm)
            plot_and_save(electronvolts_grid, self.spectrum, "Electronvolts (eV)", "Normalised Intensity", "Simulated Vibronic Emission Spectrum (Electronvolts)", "vibronic_emission_ev.png", color="deeppink", sticks=sticks_ev, xlim=xlim_ev)
        
        #Absorption
        else:
            B_si = math.pi**2 * units.SIconstants['e']**2 / (units.SIconstants['me'] * units.SIconstants['eps0'] * units.SIconstants['c'] * units.SIconstants['h']) * self.fc_options["f"] / (self.fc_options["dEH"] / units.time['s'])

            ein_coeff_abs = 2 * math.pi**2 / units.constants["c0"] * self.fc_options["f"] / self.fc_options["dEH"] * units.time['s'] / units.mass['kg']
            
            print("\nLineshapes are ignored.")
            print(f"Einstein coefficient of absorption (B) in SI: {B_si:.6e} s/kg")
            print(f"Einstein coefficient of absorption (B) in au: {ein_coeff_abs:.6e} s/kg")

            cross_sec_si_max = math.pi * units.SIconstants['e']**2 / (2 * units.SIconstants['me'] * units.SIconstants['eps0'] * units.SIconstants['c']) * self.fc_options["f"] / (self.sigma / units.time['s'] * (2 * math.pi)**0.5)
            epsilon_si_max = cross_sec_si_max * units.constants['Nl'] / math.log(10)

            cross_sec_cons_au = 2 * math.pi**2 * self.fc_options["f"] / units.constants['c0']
            epsilon_au = cross_sec_cons_au * units.constants['Nl'] / math.log(10) * 1e-3 * units.length['cm']**2

            cross_sec_cons_au_max = cross_sec_cons_au / (self.sigma * (2 * math.pi)**0.5) * units.length['cm']**2
            epsilon_au_max = cross_sec_cons_au_max * 1e-3 * units.constants['Nl'] / math.log(10)

            if self.fc_options["use_omega_omegaI0"]:
                print("\nAssuming omega/omega_I0 = 1 (sharp lineshape approximation).")
                print("The absorption cross section has uniform scaling across the spectrum.")

                print(f"\nabs_cross_sec_SI: {cross_sec_si_max * 1e20:.6e} A^2")
                print(f"max_epsilon_SI: {epsilon_si_max * 10:.6e} M^-1 cm^-1")

                print(f"\nCharacteristic maximum absorption cross section: {cross_sec_cons_au_max * 1e16:.6e} A^2")
                print(f"Characteristic maximum extinction coefficient: {epsilon_au_max:.6e} M^-1 cm^-1")
                epsilon_spectrum = epsilon_au * self.spectrum

            else:            

                print("""
Including full frequency-dependent factor (omega/omega_I0).
The absorption cross section now varies at each frequency point.
No single fixed peak formula applies.""")

                omega_ratio = self.energy_grid / self.fc_options["dEH"]
                epsilon_spectrum = epsilon_au * omega_ratio * self.spectrum

            write_table("vibronic_spectrum_data.txt", "Wavenumber(cm^-1)   Wavelength(nm)  Electronvolts(eV)   Molar Extinction Coefficients(M^-1 cm^-1)", [wavenumber_grid, wavelength_grid, electronvolts_grid, epsilon_spectrum],)
            print("\nAbsorption spectrum data saved to 'vibronic_spectrum_data.txt'")

            if plot_sticks:
                write_table("vibronic_abs_stick_data.txt", "Wavenumber(cm^-1)   Wavelength(nm)  Electronvolts(eV)   Intensity(normalised)", [stick_wavenumbers, stick_wavelengths, stick_evolts, self.stick_intensities],)
                print("Stick spectrum data saved to 'vibronic_abs_stick_data.txt'")

            stick_y = numpy.array(stick_intensities_normalised) * numpy.max(epsilon_spectrum)
            sticks_cm = (stick_wavenumbers, stick_y) if plot_sticks else None
            sticks_nm = (stick_wavelengths, stick_y) if plot_sticks else None
            sticks_ev = (stick_evolts, stick_y) if plot_sticks else None

            plot_and_save(wavenumber_grid, epsilon_spectrum, "Wavenumber (cm$^{-1}$)", "Molar Extinction Coefficient (M$^{-1}$ cm$^{-1}$)", "Simulated Vibronic Absorption Spectrum", "vibronic_spectrum.png", color="darkgreen", sticks=sticks_cm, xlim=xlim_cm,)
            plot_and_save(wavelength_grid, epsilon_spectrum, "Wavelength (nm)", "Molar Extinction Coefficient (M$^{-1}$ cm$^{-1}$)", "Simulated Vibronic Absorption Spectrum (Wavelength)", "vibronic_spectrum_nm.png", color="purple", sticks=sticks_nm, xlim=xlim_nm)
            plot_and_save(electronvolts_grid, epsilon_spectrum, "Electronvolts (eV)", "Molar Extinction Coefficient (M$^{-1}$ cm$^{-1}$)", "Simulated Vibronic Absorption Spectrum (Electronvolts)", "vibronic_spectrum_ev.png", color="brown", sticks=sticks_ev, xlim=xlim_ev)


    def run_plot(self):
        self.prep()
        self.compute_k_max()
        self.make_energy_grid()
        self.build_spectrum()
        self.plot_postrun(self.fc_options["plot_sticks"])
