# Homework 8 · Exercise 2

[![Open HW8.ipynb in Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/mj-will/fisher-example-binder/main?filepath=HW8.ipynb)

Use LIGO Hanford strain to fit the two component masses of GW150914 and compare the likelihood and Fisher contours. The badge opens `HW8.ipynb` when this version is on the public `main` branch of `mj-will/fisher-example-binder`.

**Download your completed `.ipynb` before closing Binder.** Sessions are temporary.

To export a PDF, run the cells and save your notebook in Binder. In JupyterLab, open **File**, choose **Save and Export Notebook As**, then select **PDF via LaTeX**.

To run locally, start in the repository root with Python 3.12.

```sh
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt jupyterlab
jupyter lab HW8.ipynb
```

## Hanford data and Fourier conventions

The Gravitational Wave Open Science Center (GWOSC), a service of LIGO Laboratory, the LIGO Scientific Collaboration and the Virgo Collaboration, provides the bundled 32-second, 4096-Hz Hanford strain file. The file is unchanged from the GWTC-1-confident GW150914-v3 release (DOI 10.7935/82H3-HH23), under CC BY 4.0. SHA-256 is `c5ea87beced5094b56b4d694f3a36e6d59f74d3a46e546b980a6a199d1f252a9`. [Release](https://gwosc.org/eventapi/html/GWTC-1-confident/GW150914/v3/) · [Exact HDF5 file](https://gwosc.org/eventapi/html/GWTC-1-confident/GW150914/v3/H-H1_GWOSC_4KHZ_R1-1126259447-32.hdf5) · [Data license](https://creativecommons.org/licenses/by/4.0/)

The loader and mass conversion come from the course's `HW8_start.py`. The analysis uses the same four-second Tukey-windowed strain segment, 0.25 Hz bins, and one-sided Welch PSD estimated from the first 12 seconds in two-second segments. The strain FFT includes the sample spacing and is rephased to GPS 1126259462. The templates evaluate directly on its positive-frequency bins in strain-seconds. A time offset `t` uses `exp(-2πift)`.

The inspiral band is fixed at 30–46 Hz inclusive. The minimum Schwarzschild ISCO frequency is 46.38 Hz on the grid and stays above 46 Hz with the finite-difference steps. The lower grid edge at η = 0.19 is a search limit. Contours that reach it extend beyond the sampled range at that level.

## TaylorF2 source and preliminary fit

The standalone NumPy port in `taylorf2.py` uses the nonspinning phasing coefficients and waveform convention from LALSuite tag `lalsuite-v7.26`, commit `736838ae5f43f532121c0a5627305e276d482254`. The two source headers specify GPL-2.0-or-later. Their notices are preserved in the port and the upstream GPL version 2 text is in [COPYING](COPYING). The port evaluates the 3.5PN phase and algebraically simplified leading-order stationary-phase amplitude directly on the supplied frequencies. This subset fixes the reference at formal coalescence and omits spins, tides, higher-order amplitude corrections, and LAL series allocation. The detector response and analysis band are applied in the notebook. [Upstream TaylorF2 source](https://git.ligo.org/lscsoft/lalsuite/-/blob/736838ae5f43f532121c0a5627305e276d482254/lalsimulation/lib/LALSimInspiralTaylorF2.c) · [Upstream PN coefficients](https://git.ligo.org/lscsoft/lalsuite/-/blob/736838ae5f43f532121c0a5627305e276d482254/lalsimulation/lib/LALSimInspiralPNCoefficients.c)

The waveform uses a reference distance of 400 Mpc. In the preliminary fit, the masses were fixed at the 30–350 Hz IMRPhenomD maximum from the handout's setup. The time, amplitude scale and Fourier phase were fitted to Hanford strain on 30–46 Hz and then held fixed for the exercise.

| Time offset (s) | Amplitude scale | Fourier phase rotation (rad) |
| ---: | ---: | ---: |
| 0.387507505683 | 0.398852643467 | 2.815001078781 |
