function adaptElastixTransFile(filenameIn, filenameOut, parsToAdapt, values)
%parsToAdapt = InitialTransform = InitialTransformParametersFileName
%values = C:\/....transformParameters.0.txt

if ~iscell(parsToAdapt)
    parsToAdapt = {parsToAdapt};
end
if ~iscell(values)
    values = {values};
end
% Open file and read content in string
fileIDorig  = fopen(filenameIn,'r');
stringFile = 	textscan(fileIDorig,'%s','Delimiter','','endofline','');
stringFile = stringFile{1}{1};
fclose(fileIDorig);
stringFileNew = stringFile;
for n=1:length(parsToAdapt)
    % Find location of the parameter to adapt and change
    [~,endInd] = regexp(stringFileNew,['(' parsToAdapt{n} ' ']);
    newlines = strfind(stringFileNew,char(10));
    indNewLine = find(newlines>endInd,1);
    stringFileA = stringFileNew(1:endInd-1);
    stringFileB = stringFileNew(newlines(indNewLine)-2:end); %used to be
    if isempty(str2num(values{n})) % if it's a character add double quotes
        stringFileA = [stringFileA ' "' values{n} '"'];
    else
        stringFileA = [stringFileA ' ' values{n}];
    end
    stringFileNew =[stringFileA stringFileB];
    %     stringFileNew = strrep(stringFileNew,'\','\\');
    % Save file with adapted parameter
    fileIDtemp  = fopen(filenameOut,'w');
    fprintf(fileIDtemp,'%s',stringFileNew);
    fclose(fileIDtemp);
end

% replace any possible double quotes (Linux issue)
fid  = fopen(filenameIn,'rt');
X = fread(fid);
fclose(fid);
X = char(X.');
Y = strrep(X, '""', '"');

fid2 = fopen(filenameIn,'wt');
fwrite(fid2,Y);
fclose(fid2);
end