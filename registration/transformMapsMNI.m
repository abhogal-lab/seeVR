function [] = transformMapsMNI(mapdir, opts)
global opts

elastixroot = opts.elastixdir;

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