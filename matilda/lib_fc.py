"""
version 1.0.0
author: Sayan Ghosh, Felix Plasser, Giulia Woelfle-Conway
usage: Routines for Franck-Condon simulations
"""

import math
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
        for line in open(self.hr_data_file, 'r').readlines()[1:]:
            parts = line.strip().split()
            if len(parts) < 3:
                continue
            freq_cm = float(parts[1])
            S = float(parts[2])

            freq_au = freq_cm / units.energy['rcm']
            self.freqs_au.append(freq_au)
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
                self.FC_modes += [(S, omega)]
            else:
                self.sig_modes += [(S, omega)]

        if len(self.FC_modes) == 0:
            raise error_handler.MsgError("No modes above the given S_min and w_min for FC progression")

        self.sigma = math.sqrt(sum((omega**2) * S for S, omega in self.sig_modes))

    def prep(self):
        self.read_hr_data()
        self.classify_modes()
    
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
            E_min = self.fc_options["dEH"] - (10000.0/units.energy['rcm'])      #It is 10000 cm^-1
            E_max = self.fc_options["dEH"] + (10000.0/units.energy['rcm'])
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
        self.spectrum = numpy.zeros_like(self.energy_grid)
        S_total = sum(S for S, _ in self.FC_modes)
        prefactor = math.exp(-S_total)

        def loop_over_modes(idx, E_shift, intensity):
            if idx == len(self.FC_modes):
                E_transition = self.fc_options["dEH"] + E_shift
                self.spectrum[:] += intensity * self.gaussian(self.energy_grid, E_transition)

                self.stick_energies.append(E_transition)
                self.stick_intensities.append(intensity)
                return
                
            S, omega = self.FC_modes[idx]
            for k in range(self.k_max + 1):
                new_intensity = intensity * (S**k) / math.factorial(k)
                if self.fc_options["abs_emi"] == 2:
                    new_Eshift = E_shift - k * omega
                else:
                    new_Eshift = E_shift + k * omega
                loop_over_modes(idx + 1, new_Eshift, new_intensity)

        loop_over_modes(0, 0.0, prefactor)

    def plot_postrun(self, plot_sticks: bool = True):
        """
        Prints out the results and saves data files
        """
        def write_table(filename: str, header: str, cols):
            with open(filename, "w") as fh:
                fh.write(header + "\n")
                for row in zip(*cols):
                    fh.write(" ".join(f"{float(v):15.8f}" for v in row) + "\n")
        
        def plot_and_save(x, y, xlabel, ylabel, title, filename, color, xlim, sticks=None,):
            plt.figure(figsize=(8, 5))
            plt.plot(x, y, color=color, label="Broadened spectrum")
            plt.xlim(*xlim)
            if sticks is not None:
                sx, sy = sticks
                plt.vlines(sx, 0, sy, color="blue", linewidth=1, label="Stick spectrum")
                plt.legend()
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
        stick_wavelengths = units.energy['nm'] / numpy.array(self.stick_energies)
        stick_evolts = numpy.array(self.stick_energies) * units.energy['eV']

        xlim_rcm = (float(wavenumber_grid.min()), float(wavenumber_grid.max()))
        xlim_nm = (float(wavelength_grid.min()), float(wavelength_grid.max()))
        xlim_ev = (float(electronvolts_grid.min()), float(electronvolts_grid.max()))

        #Emission
        if self.fc_options["abs_emi"] == 2:
            ein_coeff_emi = 2 * self.fc_options["dEH"]**2 * self.fc_options["f"] / units.constants['c0']**3 / units.time['s'] * self.fc_options["n_r"]**3
            lifetime = 1e9 / ein_coeff_emi

            print("\nLineshapes are ignored.")
            print(f"Einstein coefficient of spontaneous emission (A): {ein_coeff_emi:.6e} s^-1")
            print(f"Excited state lifetime: {lifetime:.6e} ns")

            if self.fc_options["use_omega_omegaI0"]:
                print("\nAssuming omega/omega_I0 = 1 (sharp lineshape approximation).")
                emission_spectrum = ein_coeff_emi * self.spectrum

            else:
                print("Including full frequency-dependent factor.")
                omega_ratio = (self.energy_grid / self.fc_options["dEH"])**3
                epsilon_spectrum = ein_coeff_emi * omega_ratio * self.spectrum

            if self.fc_options["int_type"]:
                emission_spectrum /= numpy.max(emission_spectrum)

                write_table("vibronic_spectrum_emi_normalised_data.txt", "Wavenumber(cm^-1)   Wavelength(nm)  Electronvolts(eV)   Intensity(normalised)", [wavenumber_grid, wavelength_grid, electronvolts_grid, emission_spectrum],)
                print("\nEmission spectrum data saved to 'vibronic_spectrum_emi_normalised_data.txt'")

                if plot_sticks:
                    write_table("vibronic_emi_stick_data.txt", "Wavenumber(cm^-1)  Wavelength(nm)   Electronvolts(eV)   Intensity(normalised)", [stick_wavenumbers, stick_wavelengths, stick_evolts, self.stick_intensities],)
                    print("Stick spectrum data saved to vibronic_emi_stick_data_norm.txt'")

                sticks_rcm = (stick_wavenumbers, stick_intensities_normalised) if plot_sticks else None
                sticks_nm = (stick_wavelengths, stick_intensities_normalised) if plot_sticks else None
                sticks_ev = (stick_evolts, stick_intensities_normalised) if plot_sticks else None

                plot_and_save(wavenumber_grid, emission_spectrum, "Wavenumber (cm$^{-1}$)", "Normalised Intensity", "Simulated Vibronic Emission Spectrum (Wavenumber)", "vibronic_emission_rcm_norm.png", color="darkred", sticks=sticks_rcm, xlim=xlim_rcm)
                plot_and_save(wavelength_grid, emission_spectrum, "Wavelength (nm)", "Normalised Intensity", "Simulated Vibronic Emission Spectrum (Wavelength)", "vibronic_emission_nm_norm.png", color="orange", sticks=sticks_nm, xlim=xlim_nm)
                plot_and_save(electronvolts_grid, emission_spectrum, "Electronvolts (eV)", "Normalised Intensity", "Simulated Vibronic Emission Spectrum (Electronvolts)", "vibronic_emission_ev_norm.png", color="deeppink", sticks=sticks_ev, xlim=xlim_ev)

            else:
                write_table("vibronic_spectrum_emi_rcm.txt", "Wavenumber(cm^-1) Radiative decay rate(s^-1 cm)", [wavenumber_grid, emission_spectrum/units.energy['rcm']],)
                write_table("vibronic_spectrum_emi_nm.txt", "Wavelength(nm) Radiative decay rate(s^-1 nm^-1)", [wavelength_grid, emission_spectrum/units.energy['nm']],)
                write_table("vibronic_spectrum_emi_nm.txt", "Electronvolts(eV) Radiative decay rate(s^-1 eV^-1)", [electronvolts_grid, emission_spectrum/units.energy['eV']],)
                print("\nEmission spectrum data saved to 'vibronic_spectrum_emi.txt' in cm^-1, nm and eV")

                if plot_sticks:
                    write_table("vibronic_emi_stick_data.txt", "Wavenumber(cm^-1)  Wavelength(nm)   Electronvolts(eV)   Intensity(normalised)", [stick_wavenumbers, stick_wavelengths, stick_evolts, self.stick_intensities],)
                    print("Stick spectrum data saved to vibronic_emi_stick_data.txt'")

                stick_y = numpy.array(stick_intensities_normalised) * numpy.max(emission_spectrum)
                sticks_rcm = (stick_wavenumbers, stick_y/units.energy['rcm']) if plot_sticks else None
                sticks_nm = (stick_wavelengths, stick_y/units.energy['nm']) if plot_sticks else None
                sticks_ev = (stick_evolts, stick_y/units.energy['eV']) if plot_sticks else None

                plot_and_save(wavenumber_grid, emission_spectrum/units.energy['rcm'], "Wavenumber (cm$^{-1}$)", "Radiative rate density (s$^{-1}$ cm)", "Simulated Vibronic Emission Spectrum (Wavenumber)", "vibronic_emission_rcm.png", color="darkred", sticks=sticks_rcm, xlim=xlim_rcm)
                plot_and_save(wavelength_grid, emission_spectrum/units.energy['nm'], "Wavenumber (nm)", "Radiative rate density (s$^{-1}$ nm$^{-1}$)", "Simulated Vibronic Emission Spectrum (Wavelength)", "vibronic_emission_nm.png", color="orange", sticks=sticks_nm, xlim=xlim_nm)
                plot_and_save(electronvolts_grid, emission_spectrum/units.energy['eV'], "Electronvolts (eV)", "Radiative rate density (s$^{-1}$ eV$^{-1}$)", "Simulated Vibronic Emission Spectrum (Electronvolts)", "vibronic_emission_ev.png", color="deeppink", sticks=sticks_ev, xlim=xlim_ev)

        #Absorption
        else:
            ein_coeff_abs = 2 * math.pi**2 / units.constants["c0"] * self.fc_options["f"] / self.fc_options["dEH"] / self.fc_options["n_r"] * units.time['s'] / units.mass['kg']
            
            print("\nLineshapes are ignored.")
            print(f"Einstein coefficient of absorption (B): {ein_coeff_abs:.6e} s/kg")

            cross_sec_cons_au = 2 * math.pi**2 * self.fc_options["f"] / units.constants['c0'] / self.fc_options["n_r"]
            epsilon_au = cross_sec_cons_au * units.constants['Nl'] / math.log(10) * 1e-3 * units.length['cm']**2

            cross_sec_cons_au_max = cross_sec_cons_au / (self.sigma * (2 * math.pi)**0.5) * units.length['cm']**2
            epsilon_au_max = cross_sec_cons_au_max * 1e-3 * units.constants['Nl'] / math.log(10)

            if self.fc_options["use_omega_omegaI0"]:
                print("\nAssuming omega/omega_I0 = 1 (sharp lineshape approximation).")
                print("The absorption cross section has uniform scaling across the spectrum.")

                print(f"\nCharacteristic maximum absorption cross section: {cross_sec_cons_au_max * 1e16:.6e} A^2")
                epsilon_spectrum = epsilon_au * self.spectrum

            else:            

                print("""
Including full frequency-dependent factor (omega/omega_I0).
The absorption cross section now varies at each frequency point.
No single fixed peak formula applies.""")

                omega_ratio = self.energy_grid / self.fc_options["dEH"]
                epsilon_spectrum = epsilon_au * omega_ratio * self.spectrum

            if self.fc_options["int_type"]:
                epsilon_spectrum /= numpy.max(epsilon_spectrum)
                print("\nIntensity normalized to the range [0, 1].")

                write_table("vibronic_spectrum_abs_normalized_data.txt", "Wavenumber(cm^-1)   Wavelength(nm)  Electronvolts(eV)   Normalised Intensity", [wavenumber_grid, wavelength_grid, electronvolts_grid, epsilon_spectrum])
                print("Normalized absorption spectrum data saved to 'vibronic_spectrum_abs_normalized_data.txt'")

                if plot_sticks:
                    write_table("vibronic_abs_stick_data.txt", "Wavenumber(cm^-1)   Wavelength(nm)  Electronvolts(eV)   Intensity(normalised)", [stick_wavenumbers, stick_wavelengths, stick_evolts, self.stick_intensities],)
                    print("Stick spectrum data saved to 'vibronic_abs_stick_data.txt'")

                sticks_rcm = (stick_wavenumbers, stick_intensities_normalised) if plot_sticks else None
                sticks_nm = (stick_wavelengths, stick_intensities_normalised) if plot_sticks else None
                sticks_ev = (stick_evolts, stick_intensities_normalised) if plot_sticks else None

                plot_and_save(wavenumber_grid, epsilon_spectrum, "Wavenumber (cm$^{-1}$)", "Normalised Intensity", "Simulated Vibronic Absorption Spectrum (Wavenumber)", "vibronic_absorption_rcm_norm.png", color="darkgreen", sticks=sticks_rcm, xlim=xlim_rcm)
                plot_and_save(wavelength_grid, epsilon_spectrum, "Wavelength (nm)", "Normalised Intensity", "Simulated Vibronic Absorption Spectrum (Wavelength)", "vibronic_absorption_nm_norm.png", color="purple", sticks=sticks_nm, xlim=xlim_nm)
                plot_and_save(electronvolts_grid, epsilon_spectrum, "Electronvolts (eV)", "Normalised Intensity", "Simulated Vibronic Absorption Spectrum (Electronvolts)", "vibronic_absorption_ev_norm.png", color="brown", sticks=sticks_ev, xlim=xlim_ev)

            else:
                write_table("vibronic_spectrum_abs_data.txt", "Wavenumber(cm^-1)   Wavelength(nm)  Electronvolts(eV)   Molar Extinction Coefficients(M^-1 cm^-1)", [wavenumber_grid, wavelength_grid, electronvolts_grid, epsilon_spectrum],)
                print("\nAbsorption spectrum data saved to 'vibronic_spectrum_abs_data.txt'")

                if plot_sticks:
                    write_table("vibronic_abs_stick_data.txt", "Wavenumber(cm^-1)   Wavelength(nm)  Electronvolts(eV)   Intensity(normalised)", [stick_wavenumbers, stick_wavelengths, stick_evolts, self.stick_intensities],)
                    print("Stick spectrum data saved to 'vibronic_abs_stick_data.txt'")

                stick_y = numpy.array(stick_intensities_normalised) * numpy.max(epsilon_spectrum)
                sticks_rcm = (stick_wavenumbers, stick_y) if plot_sticks else None
                sticks_nm = (stick_wavelengths, stick_y) if plot_sticks else None
                sticks_ev = (stick_evolts, stick_y) if plot_sticks else None

                plot_and_save(wavenumber_grid, epsilon_spectrum, "Wavenumber (cm$^{-1}$)", "Molar Extinction Coefficient (M$^{-1}$ cm$^{-1}$)", "Simulated Vibronic Absorption Spectrum (Wavenumber)", "vibronic_absorption_rcm.png", color="darkgreen", sticks=sticks_rcm, xlim=xlim_rcm)
                plot_and_save(wavelength_grid, epsilon_spectrum, "Wavelength (nm)", "Molar Extinction Coefficient (M$^{-1}$ cm$^{-1}$)", "Simulated Vibronic Absorption Spectrum (Wavelength)", "vibronic_absorption_nm.png", color="purple", sticks=sticks_nm, xlim=xlim_nm)
                plot_and_save(electronvolts_grid, epsilon_spectrum, "Electronvolts (eV)", "Molar Extinction Coefficient (M$^{-1}$ cm$^{-1}$)", "Simulated Vibronic Absorption Spectrum (Electronvolts)", "vibronic_absorption_ev.png", color="brown", sticks=sticks_ev, xlim=xlim_ev)

    def run_plot(self):
        self.prep()
        self.compute_k_max()
        self.make_energy_grid()
        self.build_spectrum()
        self.plot_postrun(self.fc_options["plot_sticks"])
