"""Data loading and mass conversion adapted from the course HW8_start.py."""

import h5py
import numpy as np
from scipy.signal import welch
from scipy.signal.windows import tukey


def load_spectrum(path, band_hz):
    """Return in-band frequencies, strain FFT, one-sided PSD, and df."""
    with h5py.File(path, "r") as source:
        strain = source["strain/Strain"][:]
        start = float(source["meta/GPSstart"][()])
        dt = float(source["strain/Strain"].attrs["Xspacing"])

    sample_rate = round(1 / dt)
    if start != 1126259447 or sample_rate != 4096 or len(strain) != 32 * sample_rate:
        raise ValueError("Use the 32-second, 4096-Hz GW150914-v3 HDF5 file.")

    noise = strain[: 12 * sample_rate]
    psd_f, psd = welch(noise, fs=sample_rate, nperseg=2 * sample_rate)

    first = 13 * sample_rate  # GPS 1126259460
    segment = strain[first : first + 4 * sample_rate]
    window = tukey(len(segment), alpha=0.2)
    frequencies = np.fft.rfftfreq(len(segment), d=dt)
    # Refer the Fourier phase to GPS 1126259462, the middle of the segment.
    spectrum = dt * np.fft.rfft(window * (segment - segment.mean()))
    spectrum *= np.exp(2j * np.pi * frequencies * 2)
    low, high = band_hz
    band = (frequencies >= low) & (frequencies <= high)
    in_band = frequencies[band]
    return in_band, spectrum[band], np.interp(in_band, psd_f, psd), frequencies[1]


def component_masses(chirp_mass, eta):
    """Return (m1, m2) in solar masses, with m1 >= m2. Accepts arrays."""
    total = chirp_mass * eta ** (-3 / 5)
    asymmetry = np.sqrt(1 - 4 * eta)
    return 0.5 * total * (1 + asymmetry), 0.5 * total * (1 - asymmetry)
