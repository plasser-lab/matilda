"""
version 1.0
author: Felix Plasser
Library containing unit conversion factors. Magnitudes correspond to 1 atomic unit, arranged in dictionaries.
"""

energy={'au':1.0,
        'eV':27.2113961,
        'kJ/mol':2.625500E3,
        'J':4.359744E-18,
        'kcal/mol':627.51,
        'rcm':219474.631,
        'nm':45.5633515,
        'Hz':6.579684E15,
        'K':3.157746E+05,
        'RT':3.157746E+05/298.
}
length={'au':1.0,
        'A':0.529177249,
        'm':0.529177249E-10,
        'cm':0.529177249E-8
}
time={'au':1.0,
      'fs':0.0241888432,
      's':0.0241888432E-15
}

mass={'au':1.0,
      'kg':9.10938188E-31, # electron mass
      'amu':0.000548579893 # atomic mass unit in a.u. = constants['Nl']*mass['kg']*1000
}
dipole={'D':2.54174619,
        'Cm':8.47835267E-30
}
constants={'Nl':6.02214179E23,
           'c0':137.035999, # speed of light in a.u., 1/alpha
           'h' :6.62607015E-34
}
SIconstants={'e':1.602176634E-19, # elementary charge C
           'me':9.1093837139E-31, # electron mass kg
           'c' :299792458, # speed of light in m/s
           'h' :6.62607015E-34, # planck constant J.s
           'eps0' :8.8541878188E-12  # electrical permittivity of free space F/m
}


# short cuts
def eV2nm(val):
    return energy['eV']*energy['nm']/val
