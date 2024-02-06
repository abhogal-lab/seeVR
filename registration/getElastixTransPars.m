function pars = getElastixTransPars(filenameIn,parToRead)

if ~iscell(parToRead)
    parToRead = {parToRead};
end

% Open file and read content in string
fileIDorig  = fopen(filenameIn,'r');
stringFile = 	textscan(fileIDorig,'%s','Delimiter','','endofline','');
stringFile = stringFile{1}{1};
fclose(fileIDorig);

% Find location of the parameter to adapt and change 
[~,endInd] = regexp(stringFile,['(' parToRead{1} ' ']);
newlines = strfind(stringFile,char(10));
indNewLine = find(newlines>endInd,1);
pars = str2num(stringFile(endInd:newlines(indNewLine)-3));
end