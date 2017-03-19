function filenameSysParams=setSysParamsFilename(filenameParameters)
    k=strfind(filenameParameters,'/');
    if(isempty(k))
        filenameSysParams='SysParams.mat';
    else
        k=k(length(k));
        filenameSysParams=filenameParameters(1:k);
        filenameSysParams=[filenameSysParams,'SysParams.mat'];
    end
end