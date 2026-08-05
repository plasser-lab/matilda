#!/usr/bin/python3

from __future__ import print_function

from matilda import input_options, units, lib_fc, matilda_header, error_handler
import os, sys

"""
author: Felix Plasser, Sayan Ghosh, Giulia Woelfle-Conway
usage: Spectrum using Franck-Condon progression formula based on Huang-Rhys factors.
"""

class fc_options(input_options.write_options):
    def fc_input(self):
        
        self.choose_list("Select mode", "abs_emi", [ (1, "Absorption"), (2, "Emission")], 1)
        self.read_float("Electronic adiabatic energy in Hartree", "dEH", 1.0)
        self.read_float("Oscillator strength", "f", 1.0)
        self.read_float("Refractive index of the solvent", "n_r", 1.0)
        self.read_float("Minimum vibrational frequency in cm^-1", "w_min", 200.0)
        self.read_float("Minimum Huang-Rhys factor value", "S_min", 0.1)
        self.read_yn("Show stick spectrum overlay", "plot_sticks", False)
        self.read_int("Number of points in spectrum", "npoints", 2000)
        self.read_yn("Nomralised Intensity", "int_type", True)
        self.read_yn("Assume omega/omega_I0 = 1 sharp lineshape approximation", "use_omega_omegaI0", True)

    def run_fc(self, hr_data_file):
        """
        Run Franck-Condon simulations.
        """
        print("\nRunning FC simulation...\n")
    
        hr = lib_fc.HRData(self, hr_data_file)
        hr.run_plot()

if __name__=='__main__':
    matilda_header.print_header('FCmini - Franck-Condon simulations')

    if len(sys.argv) < 2:
        print("Usage: fc_mini_v2.py <hr_data_file>")
        raise error_handler.MsgError('Specify Huang-Rhys data file')

    hr_data_file = sys.argv[1]

    fco = fc_options('fc.in')
    fco.fc_input()
    fco.run_fc(hr_data_file)
