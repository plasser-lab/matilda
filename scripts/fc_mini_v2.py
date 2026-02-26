#!/usr/bin/python3

from matilda import input_options, units, lib_fc, matilda_header

class fc_options(input_options.write_options):
    def fc_input(self):
        """
        Input for Franck-Condon simulations.
        """
        self.choose_list(
            "Select mode:",
            "spec_mode",
        [ (1, "Absorption"),
          (2, "Emission")
        ], 1)

        self.read_float("Electronic adiabatic energy (in Hartree)", "dEH", 0.)
        self.read_float("w_min(in cm^-1): ", "w_min", 200)

    def run_fc(self):
        """
        Run Franck-Condon simulations.
        """
        print("\nRunning FC simulation")

        E_0_cm = self['dEH'] * units.energy['rcm']
        print("The adiabatic energy is %f cm-1"%E_0_cm)

        # Call routines from lib_fc like this
        # lib_fc.build_spectrum(...)

if __name__=='__main__':
    matilda_header.print_header('FCmini - Franck-Condon simulations')

    fco = fc_options('fc.in')
    fco.fc_input()
    fco.run_fc()
