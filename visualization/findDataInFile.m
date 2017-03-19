function DATA=findDataInFile(DATAstring,DATAname)
	% This routine looks for the value of the variable (whose name is stored in DATAname) in the string
	% DATAstring. It is supposed that all varaibles appear as 'variable=value, ...'.
	lDATAname=length(DATAname);
	k=strfind(DATAstring,[DATAname,'=']);
	if(length(k)>=1) 
		l=strfind(DATAstring( (k(1)+lDATAname+1):length(DATAstring) ), ',' )+k(1)+lDATAname-1;
		% first ',' after the name of the DATA.
		if(isempty(l))
			DATA=Inf;
		else
			DATA=str2num(DATAstring( (k(1)+lDATAname+1):l(1) ));
		end
	else
		DATA=Inf;
	end
end