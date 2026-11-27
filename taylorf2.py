"""
/*
 *  Copyright (C) 2007 Jolien Creighton, B.S. Sathyaprakash, Thomas Cokelaer
 *  Copyright (C) 2012 Leo Singer, Evan Ochsner, Les Wade, Alex Nitz
 *  Assembled from code found in:
 *    - LALInspiralStationaryPhaseApproximation2.c
 *    - LALInspiralChooseModel.c
 *    - LALInspiralSetup.c
 *    - LALSimInspiralTaylorF2ReducedSpin.c
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with with program; see the file COPYING. If not, write to the
 *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 *  MA  02110-1301  USA
 */

/*
*  Copyright (C) 2011 Drew Keppel, 2012 Riccardo Sturani
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
*  MA  02110-1301  USA
*/
"""

# Modified 2026-11-27 for Homework 8. This is a standalone NumPy port of
# the circular, nonspinning 3.5PN phase and leading-order SPA amplitude.
# Upstream revision: 736838ae5f43f532121c0a5627305e276d482254.
# See README.md for the changes and COPYING for the upstream GPL text.

import numpy as np

SOLAR_MASS_SECONDS = 4.925490947641267e-6
SOLAR_MASS_METRES = 1476.6250380501247
PARSEC_METRES = 3.0856775814913673e16
REFERENCE_DISTANCE_METRES = 400e6 * PARSEC_METRES


def taylorf2(frequencies, chirp_mass, eta):
    """Return the positive-frequency strain spectrum in seconds at 400 Mpc.

    Masses are in solar masses. Use positive frequencies within the
    inspiral band. The face-on plus waveform has zero time shift and zero
    orbital phase with the phase reference at formal coalescence. A shift
    by t seconds multiplies this spectrum by exp(-2 pi i f t).
    """
    total_mass = chirp_mass * eta ** (-3 / 5)
    v = np.cbrt(np.pi * total_mass * SOLAR_MASS_SECONDS * frequencies)
    log_v = np.log(v)

    # v**k is the k/2 PN phasing correction to the leading term.
    p2 = 5 * (74.3 / 8.4 + 11 * eta) / 9
    p3 = -16 * np.pi
    p4 = 5 * (3058.673 / 7.056 + (5429 / 7) * eta + 617 * eta**2) / 72
    p5 = (5 / 9) * (772.9 / 8.4 - 13 * eta) * np.pi
    p6_log = -684.8 / 2.1
    p6 = (
        11583.231236531 / 4.694215680 - (640 / 3) * np.pi**2
        + p6_log * (np.euler_gamma + np.log(4))
        + eta * (-15737.765635 / 3.048192 + (225.5 / 1.2) * np.pi**2)
        + (76.055 / 1.728) * eta**2 - (127.825 / 1.296) * eta**3
    )
    p7 = np.pi * (770.96675 / 2.54016 + (378.515 / 1.512) * eta - (740.45 / 7.56) * eta**2)

    phase = (3 / (128 * eta * v**5)) * (
        1 + p2 * v**2 + p3 * v**3 + p4 * v**4
        + p5 * (1 + 3 * log_v) * v**5
        + (p6 + p6_log * log_v) * v**6 + p7 * v**7
    )
    # Algebraic simplification of the upstream Newtonian flux/energy amplitude.
    amplitude = (
        -np.sqrt(5 / 24) * np.pi**(-2 / 3)
        * SOLAR_MASS_METRES / REFERENCE_DISTANCE_METRES
        * SOLAR_MASS_SECONDS**(-1 / 6) * chirp_mass**(5 / 6)
        * frequencies**(-7 / 6)
    )
    return amplitude * np.exp(-1j * (phase - np.pi / 4))
