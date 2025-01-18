Registration Functions

!!! ALWAYS specify the path to the elasix binaries in the opts structure !!!
---> opts.elastixDir = '/seeVR/registration/elastix'

*registration is done via elastix and is guided by parameter files.
 The options contained within these files can be tweaked per use case (SEE ELASTIX MANUAL!).
 Elastix is great when you get it working but can also go wrong very easily. I have tried 
 to create optimized parameter files but they may not always work. Give thought to 
 interpolation parameter based on your applications.
 
*see functions for required inputs - generally brain extracted images should be used.
 
*functions take PATHS to files, rather than matlab variables as inputs.
 
!!!These can be used for testing things out!!!
 
***affineReg: Performs linear (affine) registration and inverse registration. Registration 
is guided by 'ParameterFileAf' and 'ParameterFileAf_rev'.
 
***nlinReg: Performs non-linear (b-spline) registration and inverse registration. Registration
is guided by 'ParameterFileBs' and 'ParameterFileBs_rev'.
 
!!!These can be used in pipelines after optimization!!!
 
funcToMNI: performs linear and non-linear registration of functional image to MNI space using
T1 intermediate. Performs inverse registration to structural and functional space. Applies 
inverse transforms to a series of masks and atlases.

funcToStruct: performs linear and non-linear registration of functional image to structural (T1)
image. Performs inverse registration to functional space.

mapsToStruct: same as above except you can supply a directory with parameter maps that will then
also be stransformed to structural space.  

structToMNI: performs linear and non-linear registration of stuctural (T1) image to MNI space. 
Performs inverse registration to structural and functional space. Applies inverse transforms
to a series of masks and atlases.

transformMapsMNI: This function uses transform files stored in the opts structure to transform
all maps in the specified folder to MNI space. You must run registration above (toMNI)
in order to use this.

transformMapsStruct: This function uses transform files stored in the opts structure to transform
all maps in the specified folder to structural space. You must run registration above (toStruct)
in order to use this.

motionCorrection: performs motion correction using affine transform. Registration 
is guided by 'ParameterFileAf_mcf'


For questions or help on how to optimize registration contact me directly a.bhogal@umcutrecht.nl