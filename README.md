[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5283595.svg)](https://doi.org/10.5281/zenodo.6532362)
[![Hits](https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fgithub.com%2Fabhogal-lab%2FseeVR&count_bg=%2379C83D&title_bg=%23555555&icon=&icon_color=%23E7E7E7&title=hits&edge_flat=false)](https://hits.seeyoufarm.com)

![seeVR_sml](https://user-images.githubusercontent.com/76702516/127998873-93f9a3d4-3a93-4ed0-acc9-7db1313c8988.png)

# seeVR
A toolbox for analyzing cerebro-vascular other related hemodynamic data:

SeeVR is an open-source toolbox MATLAB distributed under the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later option. Most of the optimization has been done using Matlab 2020+. SeeVR makes use of the parallel computing toolbox (workaround = replace parfor with for loops), the statistics toolbox and the wavelet toolbox. Wherever possible alternative implementations are avilable when these toolboxes are not available. If you run into any errors or bugs, have any feedback or would like to see new things implemented please don't hesitate to contact me at a.bhogal@umcutrecht.nl or via the contact form at https://www.seevr.nl/contact-collaboration/.

If you use any CVR functionality for your work (or any derivatives based on the seeVR toolbox) please cite the following paper:

Medullary vein architecture modulates the white matter BOLD cerebrovascular reactivity signal response to CO2: observations from high-resolution T2* weighted imaging at 7T, AA Bhogal - Neuroimage, 2021, https://doi.org/10.1016/j.neuroimage.2021.118771

and/or

Alex A. Bhogal. (2021). abhogal-lab/seeVR: Current version. Zenodo. https://doi.org/10.5281/zenodo.5283595

Functionality:

* basic functions for data handling (grabTimeseries, normTimeseries, chopTimeseries, demeanData, trAlign + WIP functions for loading of respiratory data for RespirAct only)
* nuissance regression and de-noising (filtRegressor, scrubData, genGS, remLV, smthData, denoiseData)
* convolution and fitting of HRF convolved functions (fitHRF, convHRF, convEXP) 
* spectral/hemodynamic analysis (basicCVR, glmCVRidx, fALFF, freqSpec, lagCVR)

For more details download the user manual from www.seeVR.nl or have a look at the example data and tutorials at https://www.seevr.nl/tutorials/. Some documentation still needs to be updated so check the updates section of the site regularly to see what I've fixed - by now the manual is probably out of date. 

The toolbox does include some dependencies, each provided with their own specific licences. These include:

* Wonsang You (2021). Easy Bandpass Filter (https://www.mathworks.com/matlabcentral/fileexchange/58219-easy-bandpass-filter)
* John D'Errico (2021). inpaint_nans (https://www.mathworks.com/matlabcentral/fileexchange/4551-inpaint_nans)
* Daniel Baboiu (2021). Legendre polynomials (https://www.mathworks.com/matlabcentral/fileexchange/28008-legendre-polynomials)2021.
* Matt J (2021). Extract linearly independent subset of matrix columns  (https://www.mathworks.com/matlabcentral/fileexchange/77437-extract-linearly-independent-subset-of-matrix-columns)
* Benjamin Kraus (2021). nanconv (https://www.mathworks.com/matlabcentral/fileexchange/41961-nanconv)
* John D'Errico (2021). SLM - Shape Language Modeling (https://www.mathworks.com/matlabcentral/fileexchange/24443-slm-shape-language-modeling)
* Oliver Woodford (2021). SC - powerful image rendering (https://github.com/ojwoodford/sc)
* Stephen Cobeldick (2021). ColorBrewer: Attractive and Distinctive Colormaps (https://github.com/DrosteEffect/BrewerMap)
* Ander Biguri (2021). Perceptually uniform colormaps (https://www.mathworks.com/matlabcentral/fileexchange/51986-perceptually-uniform-colormaps)
* Jimmy Shen (2021). Tools for NIfTI and ANALYZE image (https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image), MATLAB Central File Exchange. Retrieved August 12, 2021.
