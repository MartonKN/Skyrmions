function generateRunFile(filenameStub,filenameFilenames,filenameTmpSysParams)
% Set source folder on the cluster
sourceFolder='/home/kanasznm/Skyrmions/ImaginaryTime/currentRuns';

filenameRunFile=[filenameStub,'__RunFile'];
load(filenameTmpSysParams);
runfileID=fopen(filenameRunFile,'w');
fprintf(runfileID,'#!/bin/bash\n');
fprintf(runfileID,'#PBS -N job__%s\n',filenameStub);

% Set number of processors
if(strcmp(SysParams__nthreads,'USE_MAX_NUM_PROCS'))
    fprintf(runfileID,'#PBS -l nodes=1:ppn=%d\n',8);
else
    fprintf(runfileID,'#PBS -l nodes=1:ppn=%d\n',...
            round(str2double(SysParams__nthreads)));
end
fprintf(runfileID,'#PBS -l walltime=120:00:00\n');

% Set required memory
memory=SysParams__Mx*SysParams__My*SysParams__Mz*8*(3*4+6)...
      +SysParams__AbsPsi_Steps*SysParams__AbsF_Steps...
      *SysParams__SqrtMinusDmu_Steps*8*7; % in bytes
memory=memory*1.1+300*1024^2; % some extra memory
memory=round(memory/1024^2);
fprintf(runfileID,'#PBS -l mem=%dmb\n',memory);

% Set required filespace
diskSpace=SysParams__Mx*SysParams__My*SysParams__Mz*...
      ((length(SysParams__saving_times)+1)*9 + 6)*8; % in bytes
diskSpace=diskSpace*1.1+150*1024^2; % just some extra space to be sure
diskSpace=round(diskSpace/1024^2); % in MegaBytes
fprintf(runfileID,'#PBS -l file=%dM\n',diskSpace);

fprintf(runfileID,'#PBS -j oe\n');
fprintf(runfileID,'#PBS -o %s/%s.out\n',sourceFolder,filenameStub);
fprintf(runfileID,'\n');
fprintf(runfileID,'\n################\n');
fprintf(runfileID,['SRC_FOLDER="',sourceFolder,'"\n']);
fprintf(runfileID,'EXE="imagTimeGPESolver"\n');
fprintf(runfileID,'################\n\n');
fprintf(runfileID,'# mail is sent to you when the job starts and when it terminates or aborts\n');
fprintf(runfileID,'#PBS -m abe -M kanasz-nagy.marton@wigner.bme.hu\n');
fprintf(runfileID,'\n');
fprintf(runfileID,'cd $PBS_O_WORKDIR\n');
fprintf(runfileID,'export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/opt/intel10/lib/:/opt/mkl-10.0/lib/em64t\n');
fprintf(runfileID,'\n');
fprintf(runfileID,'# PREPARATIONS FOR RUNNING:\n');
fprintf(runfileID,'WORKDIR="WDIR__%s"\n',filenameStub);
fprintf(runfileID,'mkdir $WORKDIR\n');
fprintf(runfileID,'cd $WORKDIR\n');
fprintf(runfileID,'# copy input and exe files\n');
fprintf(runfileID,'cp $SRC_FOLDER/%s* ./\n',filenameStub);
fprintf(runfileID,'cp $SRC_FOLDER/$EXE ./\n');
fprintf(runfileID,'# Put run details to stdout:\n');
fprintf(runfileID,'# Print system hostname\n');
fprintf(runfileID,'hostname\n');
fprintf(runfileID,'# Print path to the current directory\n');
fprintf(runfileID,'pwd\n');
fprintf(runfileID,'# ... and list its content\n');
fprintf(runfileID,'ls\n');
fprintf(runfileID,'\n');
fprintf(runfileID,'# NECESSARY SCRIPT FOR COMPILATION ON SHELDON\n');
fprintf(runfileID,'module add mkl-10.0\n');
fprintf(runfileID,'# RUN PROGRAM\n');
fprintf(runfileID,['./$EXE ',filenameFilenames,'\n']);
fprintf(runfileID,'\n');
fprintf(runfileID,'# COPY OUTPUT FILES\n');
fprintf(runfileID,'cp ./* $SRC_FOLDER/\n');
fprintf(runfileID,'rm -f *\n');
fprintf(runfileID,'\n');
fprintf(runfileID,'# EXIT\n');
fprintf(runfileID,'exit 0\n');

% Finally put a string to the screen with which I can copy these files to
% the cluster
outStr='scp ';
outStr=[outStr,pwd,'/'];
outStr=[outStr,filenameStub,'* '];
outStr=[outStr,'kanasznm@sheldon.physik.fu-berlin.de:',sourceFolder,'/'];
disp(outStr);
disp('ssh kanasznm@sheldon.physik.fu-berlin.de');
disp(['qsub ',sourceFolder,'/',filenameRunFile]);
end