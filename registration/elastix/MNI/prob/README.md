# CBV BRAIN ATLAS 

This CBV atlas was generated using the CBV maps compiled by the multicenter and international retrospective clinical study NCT03439332. A total of 134 cases from this dataset with MRI perfusion studies were included. All patients were diagnosed with Glioblastoma grade IV WHO with histopathological confirmation and followed Stupp standard treatment. A summary of the MRI acquisition protocol used by each center was presented in [1]. Quantification of rCBV and relative cerebral blood flow was performed employing standard techniques proposed in the literature [2], and processed automatically using the ONCOhabitats services [3]. Gamma-variate curve-fitting and Boxerman leakage correction technique [4] were used to correct for T2 and T1 leakage-effects.
A CBV atlas was generated using the CBV maps of the patients included in the study. These maps were non-linearly co-registered with ANTS software (SyN algorithm) [5] using cross correlation as a similarity metric. Finally, the atlas was calculated using the median CBV value at each voxel, excluding the areas affected by the tumor lesion for each patient.

![alt text](https://github.com/eliesfuster/CBV_ATLAS/blob/main/Preview.png)
Figure: On top, T1w template from MNI152 atlas. On the bottom, the CBV atlas generated from the NCT03439332 study cohort.

1. Álvarez-Torres, M. D. M. et al. Robust association between vascular habitats and patient prognosis in glioblastoma: An international multicenter study. J Magn Reson Imaging (2019) doi:10.1002/jmri.26958.
2. Boxerman, J. L. et al. Consensus recommendations for a dynamic susceptibility contrast MRI protocol for use in high-grade gliomas. Neuro Oncol 22, 1262–1275 (2020).
3. Juan-Albarracín, J., Fuster-Garcia, E., García-Ferrando, G. A. & García-Gómez, J. M. ONCOhabitats: A system for glioblastoma heterogeneity assessment through MRI. International Journal of Medical Informatics 128, 53–61 (2019).
4. Boxerman, J. L., Schmainda, K. M. & Weisskoff, R. M. Relative cerebral blood volume maps corrected for contrast agent extravasation significantly correlate with glioma tumor grade, whereas uncorrected maps do not. AJNR Am J Neuroradiol 27, 859–867 (2006).
5. Avants, B. B. et al. A reproducible evaluation of ANTs similarity metric performance in brain image registration. Neuroimage 54, 2033–2044 (2011).
