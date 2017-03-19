function [Psi1,Psi2,Psi3,n11,n12,n13,n22,n23,n33,nTOF11,nTOF12,nTOF13,nTOF22,nTOF23,nTOF33,X,Y,Z]=...
          getPsisEllipsoid(filename_Psis,filename_SysParams,rx,ry,rz,...
                           N_points,Hedgehogize)
   % Returns the values of the order parameters on the
   % ellipsoid around the origin, with the equitoral radii rx and ry 
   % and the polar radius rz.
   % X,Y,Z: contain the x, y and z coordinates of the points of the 
   % ellipsoid. (Arrays of size [N_points,N_points]).
   % The order parameters were stored by the GPE solver in 'filename'.
   % Note: before the call of this function, get_SysParams() needs to be
   % called, in order that it saves the system parameters in the file
   % filename_SysParams.
   %
   % Hedgehogize
   % 0:    nothing happens
   % 1:    Psi -> |Psi|*e_r, where e_r is the radial eigenvector
   % 2:    Psi -> Psi-|Psi|*e_r
   % This, not too often used, option helps me determine how accurately
   % is our configuration a hedgehog. This is very important to know when
   % one calculates the excitation spectrum, since the program calculating
   % it assumes completely accurate hedgehog.

   load(filename_SysParams);
   i=sqrt(-1);
   
   % Setting the points of the ellipsoid
   rx=abs(rx); ry=abs(ry); rz=abs(rz);
   N_points=max([N_points,4]);
   if (rx>SysParams__Mx*SysParams__ax) 
       rx=SysParams__Mx*SysParams__ax;
       disp('Error in function get_Psis_ellipsoid:');
       disp('Ellipsoid reaches out of the simulated volume.');
       disp(['Value of rx has been modified to ',num2str(rx)]);
   end
   if (ry>SysParams__My*SysParams__ay) 
       ry=SysParams__My*SysParams__ay;
       disp('Error in function get_Psis_ellipsoid:');
       disp('Ellipsoid reaches out of the simulated volume.');
       disp(['Value of ry has been modified to ',num2str(ry)]);
   end
   if (rz>SysParams__Mz*SysParams__az) 
       rz=SysParams__Mz*SysParams__az;
       disp('Error in function get_Psis_ellipsoid:');
       disp('Ellipsoid reaches out of the simulated volume.');
       disp(['Value of rz has been modified to ',num2str(rz)]);
   end
   [X,Y,Z]=ellipsoid(0,0,0,rx,ry,rz,N_points-1);
   
   % Coordinates of the surface points in the grid 
   MX=round(X/SysParams__ax + 0.5*(SysParams__Mx+1));
   MY=round(Y/SysParams__ay + 0.5*(SysParams__My+1));
   MZ=round(Z/SysParams__az + 0.5*(SysParams__Mz+1));
   for j1=1:N_points
       for j2=1:N_points
           MX(j1,j2)=max([min([MX(j1,j2),SysParams__Mx]),1]);
           MY(j1,j2)=max([min([MY(j1,j2),SysParams__My]),1]);
           MZ(j1,j2)=max([min([MZ(j1,j2),SysParams__Mz]),1]);
       end
   end
   
   % Fill Psi-s and n-s
   Psi1=zeros(N_points,N_points);
   Psi2=zeros(N_points,N_points);
   Psi3=zeros(N_points,N_points);
   n11=zeros(N_points,N_points);
   n12=zeros(N_points,N_points);
   n13=zeros(N_points,N_points);
   n22=zeros(N_points,N_points);
   n23=zeros(N_points,N_points);
   n33=zeros(N_points,N_points);
   nTOF11=zeros(N_points,N_points);
   nTOF12=zeros(N_points,N_points);
   nTOF13=zeros(N_points,N_points);
   nTOF22=zeros(N_points,N_points);
   nTOF23=zeros(N_points,N_points);
   nTOF33=zeros(N_points,N_points);
   tmp_M=SysParams__Mx*SysParams__My*SysParams__Mz;
   fileID=fopen(regexprep(filename_Psis,' ',''));
   for j1=1:N_points
       for j2=1:N_points
           fseek(fileID,8*(MZ(j1,j2) + SysParams__Mz*((MY(j1,j2)-1)+SysParams__My*(MX(j1,j2)-1))),'bof');
           tmp_array=fread(fileID,24,'double=>double',8*(tmp_M-1)); 
           Psi1(j1,j2)=tmp_array(1)+1i*tmp_array(2);
           Psi2(j1,j2)=tmp_array(3)+1i*tmp_array(4);
           Psi3(j1,j2)=tmp_array(5)+1i*tmp_array(6);
           n11(j1,j2)=tmp_array(7);
           n12(j1,j2)=tmp_array(8)+1i*tmp_array(9);
           n13(j1,j2)=tmp_array(10)+1i*tmp_array(11);
           n22(j1,j2)=tmp_array(12);
           n23(j1,j2)=tmp_array(13)+1i*tmp_array(14);
           n33(j1,j2)=tmp_array(15);
           nTOF11(j1,j2)=tmp_array(16);
           nTOF12(j1,j2)=tmp_array(17)+1i*tmp_array(18);
           nTOF13(j1,j2)=tmp_array(19)+1i*tmp_array(20);
           nTOF22(j1,j2)=tmp_array(21);
           nTOF23(j1,j2)=tmp_array(22)+1i*tmp_array(23);
           nTOF33(j1,j2)=tmp_array(24);
       end
   end
   fclose(fileID);
   clear MX MY MZ;
   
   % Hedgehogize
   if Hedgehogize==1
       R=sqrt(X.^2+Y.^2+Z.^2);
       R=(R>0).*R+(R==0);
       PsiNorm=sqrt(abs(Psi1).^2+abs(Psi2).^2+abs(Psi3).^2);
       Psi1=PsiNorm.*X./R;
       Psi2=PsiNorm.*Y./R;
       Psi3=PsiNorm.*Z./R;
   elseif Hedgehogize==2
       R=sqrt(X.^2+Y.^2+Z.^2);
       R=(R>0).*R+(R==0);
       PsiNorm=sqrt(abs(Psi1).^2+abs(Psi2).^2+abs(Psi3).^2);
       Psi1=Psi1-PsiNorm.*X./R;
       Psi2=Psi2-PsiNorm.*Y./R;
       Psi3=Psi3-PsiNorm.*Z./R;       
   end
end