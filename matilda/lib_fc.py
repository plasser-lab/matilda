"""
Routines for Franck-Condon simulations.
"""

def build_spectrum(modes, k_max, E_0_cm, energy_grid, sigma):
    """
    Main routine for building the spectrum.
    """
    spectrum = numpy.zeros_like(energy_grid)
    S_total = sum(S for S, _ in modes)
    prefactor = math.exp(-S_total)

    # Recursive loop
    def loop_over_modes(idx, E_shift, intensity):
        if idx == len(modes):
            # All modes handled -> place peak
            E_transition = E_0_cm + E_shift
            spectrum[:] += intensity * gaussian(energy_grid, E_transition, sigma)
            return
        S, omega = modes[idx]
        for k in range(k_max + 1):
            new_intensity = intensity * (S**k) / math.factorial(k)
            if mode_type == "emi":
                new_Eshift = E_shift - k * omega
            else:
                new_Eshift = E_shift + k * omega
            loop_over_modes(idx + 1, new_Eshift, new_intensity)

    # Start recursion
    loop_over_modes(0, 0.0, prefactor)
    return spectrum
