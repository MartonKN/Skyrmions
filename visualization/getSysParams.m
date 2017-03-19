function getSysParams(filename_Parameters,filename_SysParams)
	fileID=fopen(regexprep(filename_Parameters,' ',''),'r');
%    strcmp(regexprep(filename_Parameters,' ',''),'/home/marton/Asztal/PhD/Skyrmionok/2_TestRuns/1_Testing/WDIR__diffusionOnly_tmax=100_M=64__2011_12_16__13_13/diffusionOnly_tmax=100_M=64__2011_12_16__13_13__SystemParameters'),
%    regexprep(filename_Parameters,' ','')
%    '/home/marton/Asztal/PhD/Skyrmionok/2_TestRuns/1_Testing/WDIR__diffusionOnly_tmax=100_M=64__2011_12_16__13_13/diffusionOnly_tmax=100_M=64__2011_12_16__13_13__SystemParameters',
%    fileID,
	DATAstring=fscanf(fileID,'%c',inf);
	SysParams__Mx=findDataInFile(DATAstring,'Mx');
	SysParams__My=findDataInFile(DATAstring,'My');
	SysParams__Mz=findDataInFile(DATAstring,'Mz');
	SysParams__ax=findDataInFile(DATAstring,'ax');
	SysParams__ay=findDataInFile(DATAstring,'ay');
	SysParams__az=findDataInFile(DATAstring,'az');
	SysParams__tmax=findDataInFile(DATAstring,'tmax');
	SysParams__dt=findDataInFile(DATAstring,'dt');
	SysParams__T=findDataInFile(DATAstring,'T');
	SysParams__V0=findDataInFile(DATAstring,'V');
	SysParams__V2=findDataInFile(DATAstring,'V2');
	SysParams__J=findDataInFile(DATAstring,'J');
	SysParams__mu_0=findDataInFile(DATAstring,'mu_0');
    SysParams__trapx=findDataInFile(DATAstring,'trapx');
    SysParams__trapy=findDataInFile(DATAstring,'trapy');
    SysParams__trapz=findDataInFile(DATAstring,'trapz');
    DATAname='nthreads';
    lDATAname=length(DATAname);
	k=strfind(DATAstring,[DATAname,'=']);
    l=strfind(DATAstring( (k(1)+lDATAname+1):length(DATAstring) ), ',' )+k(1)+lDATAname-1;
    SysParams__nthreads=DATAstring((k(1)+lDATAname+1):l(1));
    if(~strcmp(SysParams__nthreads,'USE_MAX_NUM_PROCS') && ...
       isempty(str2num(SysParams__nthreads)))
        SysParams__nthreads='2';
        disp(['The value of SysParams__nthreads has been modified to ',...
              SysParams__nthreads,' in function getSysParams.']);
    end
    
	j=0;
	tmp=findDataInFile(DATAstring,['saving_times[',num2str(0),']']);
	while tmp~=Inf,
		j=j+1;
		SysParams__saving_times(j)=tmp;
		tmp=findDataInFile(DATAstring,['saving_times[',num2str(j),']']);
	end
	SysParams__dispform_saving_times(1)= ...
        findDataInFile(DATAstring,'dispform_saving_times[0]');
	SysParams__dispform_saving_times(2)= ...
        findDataInFile(DATAstring,'dispform_saving_times[1]');
	fclose(fileID);
	
	j=strfind(filename_Parameters,'SystemParameters');
	filename_base=filename_Parameters(1:(j-1));
	t=0; j=1; j1_previous=-1; tdisp_previous=-1;
	format_str=['%',num2str(SysParams__dispform_saving_times(1)),'.', ...
                num2str(SysParams__dispform_saving_times(2)),'f'];
	while t<=SysParams__tmax
		for j1=1:length(SysParams__saving_times)
			if (SysParams__saving_times(j1)>=t-SysParams__dt/2 && ... 
                SysParams__saving_times(j1)<t+SysParams__dt/2 && ...
                ~(j1==j1_previous) && ~(tdisp_previous==t))
                if(t<10)
    				SysParams__filenames(j,:)= ...
                        [filename_base,'OrderParameters__',num2str(t,format_str),'  '];
                    SysParams__saving_times_in_filenames(j,:)=[num2str(t,format_str),'  '];
                elseif(t<100)
    				SysParams__filenames(j,:)= ...
                        [filename_base,'OrderParameters__',num2str(t,format_str),' '];
                    SysParams__saving_times_in_filenames(j,:)=[num2str(t,format_str),' '];
                else
    				SysParams__filenames(j,:)= ...
                        [filename_base,'OrderParameters__',num2str(t,format_str)];
                    SysParams__saving_times_in_filenames(j,:)=num2str(t,format_str);
                end
				j1_previous=j1;
				tdisp_previous=t;
				j=j+1;
			end
		end
		t=t+SysParams__dt;
	end
	
	save(filename_SysParams,...
        'SysParams__Mx','SysParams__My','SysParams__Mz','SysParams__ax',...
        'SysParams__ay','SysParams__az','SysParams__tmax','SysParams__dt',... 
        'SysParams__T','SysParams__V0','SysParams__V2',...
	    'SysParams__J','SysParams__mu_0',...
        'SysParams__trapx','SysParams__trapy','SysParams__trapz',...
        'SysParams__saving_times',... 
        'SysParams__dispform_saving_times','SysParams__filenames',...
        'SysParams__saving_times_in_filenames','SysParams__nthreads');
end
