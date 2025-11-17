[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
![MATLAB](https://img.shields.io/badge/MATLAB-R2020+-orange.svg)
[![GitHub release](https://img.shields.io/github/v/release/abhogal-lab/seeVR)](https://github.com/abhogal-lab/seeVR/releases)
[![Hits](https://hits.sh/github.com/abhogal-lab/seeVR.svg?style=flat-square&color=blue&label=views)](https://hits.sh/github.com/abhogal-lab/seeVR/)

# seeVR — cerebrovascular MRI analysis toolbox

**seeVR** is an open-source MATLAB toolbox for analyzing cerebrovascular imaging data — especially CVR (cerebrovascular reactivity) based on BOLD/fMRI. It ships with ready-to-use functions for CVR mapping, lag analysis, HRF modeling, denoising, spectral metrics (e.g., fALFF), plus a 3-D/4-D viewer (**seeVieweR**) for high-quality visualization and export.

---

<img src="docs/logo.png" alt="Logo" width="300"/>

---

## Why use seeVR?

- **Purpose-built for CVR & hemodynamics**: pipelines and helpers aimed at hypercapnia/breath-hold/resting-state use cases - goto www.seeVR.nl.  
- **From data to maps**: CVR, lag-corrected CVR, HRF/dispersion/onset metrics, spectral features (fALFF), plus utilities (regression, filtering, masks).  
- **Batteries included**: seeVieweR 3-D/4-D viewer for overlays, masks, erosion, slice controls, and movie/snapshot export.  

---

## Features

### Analysis functions
- **CVR mapping**: basic/GLM CVR indices, lag-aware CVR. (`basicCVR`, `gasslessCVRx`, `lagCVR`) 

<img src="docs/CVR.png" alt="Logo" width="800"/>

- **HRF modeling**: convolve and fit (onset/dispersion) via HRF/EXP kernels. (`fitTau`, `convHRF`, `fitHRF`)  

<img src="docs/dispersion.png" alt="Logo" width="800"/>

- **Aligment of physiological traces**: manual or automatic. (`autoAlign`, `trAlign`)  
- **Spectral metrics**: fALFF and frequency spectra. (`fALFF`, `freqSpec`)  
- **Preprocess & cleanup**: nuisance regression, scrubbing, band-pass, normalization, timeseries utilities. (`filtRegressor`, `scrubData`, `genGS`, `remLV`, `smthData`, `denoiseData`, `grabTimeseries`, `normTimeseries`, `chopTimeseries`, `demeanData`, `trAlign`)  
- **I/O & helpers**: NIfTI handling (via built-in functions or Jimmy Shen’s toolbox), plotting, registration helpers, and more.  

### Visualization (seeVieweR)
- **3-D/4-D volume viewer** with overlays, color/alpha control, slice masks (X/Y/Z), vessel erosion (skeleton protection), colormap flipping, lighting/denoising toggles, background, and export (PNG/SVG/GIF/MP4).  

![LOGO](docs/volumetric.png)

### Tutorials & docs
- Walk-throughs for **CO₂ CVR** and **breath-hold/resting-state CVR**, with template scripts.

Tutorials and test data available at the [seeVR toolbox website](https://www.seevr.nl/tutorials/)

---

### Bug Reporting
- contact me if you need some advice on how to put together a pipeline
- functions have been developped on a Windows machine, they probably work in Linux environment, but no guarantees for MAC. That said, I will work through errors in any case.

---

## Repository layout (high level)

- `main_functions/`, `functions/`, `helper/`, `plotting/`, `read_write/`, `registration/` — analysis code, utilities & helpers.  
- `seeVieweR/` — [Visualization tool](https://github.com/abhogal-lab/seeVieweR) (3-D/4-D).  
- `gas_import/`, `resources/` — domain-specific helpers/resources.  

---

## Installation

### Option A — clone for the toolbox (analysis code)
```bash
git clone https://github.com/abhogal-lab/seeVR.git

- Or just download and add to MATLAB path 

