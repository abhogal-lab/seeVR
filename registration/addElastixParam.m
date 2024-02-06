function  addElastixParam(filenameIn, filenameOut, paramString)

fileID = fopen(filenameIn, 'a');
fwrite(fileID, paramString);
disp(['appending: ',paramString])
fclose(fileID)

end