function [Psi1,Psi2,Psi3,n11,n12,n13,n22,n23,n33,nTOF11,nTOF12,nTOF13,nTOF22,nTOF23,nTOF33,x,y,z]=...
         getPsisInterspaced(filename_Psis,filename_SysParams,...
                            intersp_x,intersp_y,intersp_z,Hedgehogize)
   % Returns order parameter values Psi1 and Psi2 stored by the GPE solver 
   % in the file 'filename_Psis'. Not all values are returned, only 
   % each 'intersp_x'th in the x direction, each 'intersp_y'th in the
   % y direction etc. The 1D arrays x, y, z store the coordinates 
   % corresponding to these points.
   % Note: before the call of this function, get_SysParams() needs to be
   % called, in order that it saves the system parameters in the file
   % 'file_SysParams'.
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
   
   % intersp_... need to be integers within 0 and SysParams__M...
   intersp_x=max([min([floor(intersp_x),SysParams__Mx]),0]);
   intersp_y=max([min([floor(intersp_y),SysParams__My]),0]);
   intersp_z=max([min([floor(intersp_z),SysParams__Mz]),0]);
   
   % Fill Psi1 and Psi2
   fileID=fopen(regexprep(filename_Psis,' ',''));
   Psi1=zeros(floor(SysParams__Mx/intersp_x), ...
              floor(SysParams__My/intersp_y), ...
              floor(SysParams__Mz/intersp_z));
   Psi2=zeros(floor(SysParams__Mx/intersp_x), ...
              floor(SysParams__My/intersp_y), ...
              floor(SysParams__Mz/intersp_z));
   Psi3=zeros(floor(SysParams__Mx/intersp_x), ...
              floor(SysParams__My/intersp_y), ...
              floor(SysParams__Mz/intersp_z));
   n11=Psi1;
   n12=Psi1;
   n13=Psi1;
   n22=Psi1;
   n23=Psi1;
   n33=Psi1;
   nTOF11=Psi1;
   nTOF12=Psi1;
   nTOF13=Psi1;
   nTOF22=Psi1;
   nTOF23=Psi1;
   nTOF33=Psi1;
   tmp_Mx=Psi1;
   
   
   % Fill arrays
   for j=1:24
       for jx=1:(SysParams__Mx/intersp_x)
          for jy=1:(SysParams__My/intersp_y)
             fseek(fileID,8*SysParams__Mz*((jy*intersp_y-1)+(jx*intersp_x-1+(j-1)*SysParams__Mx)*SysParams__My),'bof');
             tmp_array(1,1,:) = fread(fileID,floor(SysParams__Mz/intersp_z),...
                                      'double=>double',8*(intersp_z-1)); 
             tmp_Mx(jx,jy,:)=tmp_array;
          end
       end
       switch j
           case 1
               Psi1=tmp_Mx;
           case 2
               Psi1=Psi1+1i*tmp_Mx;
           case 3
               Psi2=tmp_Mx;
           case 4
               Psi2=Psi2+1i*tmp_Mx;
           case 5
               Psi3=tmp_Mx;
           case 6
               Psi3=Psi3+1i*tmp_Mx;
           case 7
               n11=tmp_Mx;
           case 8
               n12=tmp_Mx;
           case 9
               n12=n12+1i*tmp_Mx;
           case 10
               n13=tmp_Mx;
           case 11
               n13=n13+1i*tmp_Mx;
           case 12
               n22=tmp_Mx;
           case 13
               n23=tmp_Mx;
           case 14
               n23=n23+1i*tmp_Mx;
           case 15
               n33=tmp_Mx;
           case 16
               nTOF11=tmp_Mx;
           case 17
               nTOF12=tmp_Mx;
           case 18
               nTOF12=nTOF12+1i*tmp_Mx;
           case 19
               nTOF13=tmp_Mx;
           case 20
               nTOF13=nTOF13+1i*tmp_Mx;
           case 21
               nTOF22=tmp_Mx;
           case 22
               nTOF23=tmp_Mx;
           case 23
               nTOF23=nTOF23+1i*tmp_Mx;
           case 24
               nTOF33=tmp_Mx;
       end
   end
   fclose(fileID); 
   clear tmp_array;
   
   % Fill up x, y and z
   x=((1:intersp_x:SysParams__Mx)-0.5*(1+SysParams__Mx))*SysParams__ax;
   y=((1:intersp_y:SysParams__My)-0.5*(1+SysParams__My))*SysParams__ay;
   z=((1:intersp_z:SysParams__Mz)-0.5*(1+SysParams__Mz))*SysParams__az;
   
   % Hedgehogize
   if Hedgehogize==1 || Hedgehogize==2
       X=ones(floor(SysParams__Mx/intersp_x), ...
              floor(SysParams__My/intersp_y), ...
              floor(SysParams__Mz/intersp_z));
       Y=ones(floor(SysParams__Mx/intersp_x), ...
              floor(SysParams__My/intersp_y), ...
              floor(SysParams__Mz/intersp_z));
       Z=ones(floor(SysParams__Mx/intersp_x), ...
              floor(SysParams__My/intersp_y), ...
              floor(SysParams__Mz/intersp_z));
       for i1=1:floor(SysParams__Mx/intersp_x)
           X(i1,:,:)=X(i1,:,:)*x(i1);
       end
       for i1=1:floor(SysParams__My/intersp_y)
           Y(:,i1,:)=Y(:,i1,:)*y(i1);
       end
       for i1=1:floor(SysParams__Mz/intersp_z)
           Z(:,:,i1)=Z(:,:,i1)*z(i1);
       end
       R=sqrt(X.^2+Y.^2+Z.^2);
       R=(R>0).*R+(R==0);
       PsiNorm=sqrt(abs(Psi1).^2+abs(Psi2).^2+abs(Psi3).^2);
   end
   if Hedgehogize==1
       Psi1=PsiNorm.*X./R;
       Psi2=PsiNorm.*Y./R;
       Psi3=PsiNorm.*Z./R;
   elseif Hedgehogize==2
       Psi1=Psi1-PsiNorm.*X./R;
       Psi2=Psi2-PsiNorm.*Y./R;
       Psi3=Psi3-PsiNorm.*Z./R;       
   end
   clear X Y Z R PsiNorm;
end