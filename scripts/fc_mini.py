#!/usr/bin/python3

from __future__ import print_function

from matilda import input_options, units, lib_fc, matilda_header, error_handler
import os, sys

"""
version 1.0.0
author: Felix Plasser
usage: Spectrum using Franck-Condon progression formula based on Huang-Rhys factors.
"""

class fc_options(input_options.write_options):
    def fc_input(self):
        
        self.choose_list("Select mode:", "abs_emi", [ (1, "Absorption"), (2, "Emission")], 1)
       
        #self.abs_emi = self["abs_emi"]

        self.read_float("Electronic adiabatic energy (in Hartree): ", "dEH", 1.0)
        self.read_float("Oscillator strength (dimensionless): ", "f", 1.0)
        self.read_float("w_min(in cm^-1): ", "w_min", 200.0)
        self.read_float("S_min: ", "S_min", 0.1)
        self.read_yn("Show stick spectrum overlay", "plot_sticks", False)

        self.read_float("Number of points in spectrum", "npoints", 2000)
    
        if self["abs_emi"] == 1:
            self.read_yn("Assume omega/omega_I0 = 1 sharp lineshape approximation for the absorption cross-section?", "use_omega_omegaI0", True)   

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
        raise error_handler.MsgError('Specify Huang Rhys data file')

    hr_data_file = sys.argv[1]

    fco = fc_options('fc.in')
    fco.fc_input()
    fco.run_fc(hr_data_file)
