function [] = transformMapsMNI(mapdir, opts)
global opts

elastixroot = opts.elastixdir;

%if bspline parameter file is not provided, check to see if we can find it
if ~isfield(opts,'bsplineTxParamFile')
    disp('...checking for parameter file')
    if exist(fullfile(opts.savedir,'structToMNI/mTransformParameters.1.txt'))
    opts.bsplineTxParamFile = fullfile(opts.savedir,'structToMNI/mTransformParameters.1.txt')
    else
        error('run registration to generate transformation files')
    end
end
        
%setup OS-dependent paths
if ispc
    elastixrootOS = fullfile(elastixroot,'windows');
elseif ismac
    elastixrootOS = fullfile(elastixroot,'mac','bin');
else
    elastixrootOS = fullfile(elastixroot,'linux','bin');
end
disp(['path to parameter maps: ',mapdir])
dirinfo = dir(mapdir);

opts.regOutput = fullfile(mapdir,'regOutput'); mkdir(opts.regOutput);
disp(['transformed parameter maps will be saved in: ', opts.regOutput])

%apply Tx
for kk=1:size(dirinfo,1)
    if dirinfo(kk).isdir
    else
        TxImg = fullfile(dirinfo(kk).folder,dirinfo(kk).name);
        [~] = transformixReg(TxImg, opts.bsplineTxParamFile, opts.regOutput, elastixrootOS);

     end
end