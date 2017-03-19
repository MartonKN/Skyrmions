function [Psi1,Psi2,Psi3,n11,n12,n13,n22,n23,n33,nTOF11,nTOF12,nTOF13,nTOF22,nTOF23,nTOF33,x,z]=...
          getPsisFixedY(filename_Psis,filename_SysParams,y_value,Hedgehogize)
   % Returns the values of the order parameters in the plane
   % with fixed 'y_value' distance from the origin. 
   % x and z are 1D arrays, containing the coordinates of the mesh.
   % The order parameters were stored by the GPE solver in 'filename'.
   % Note: before the call of this function, get_SysParams() needs to be
   % called, in order that it saves the system parameters in the file
   % 'filename_SysParams'.
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
   my_value=min([max([round(y_value/SysParams__ay+0.5*(1+SysParams__My)),1]), ...
                 SysParams__My]);
   i=sqrt(-1);
   fileID=fopen(regexprep(filename_Psis,' ',''));
   % Fill RePsi1
   Psi1=zeros(SysParams__Mx,SysParams__Mz);
   Psi2=zeros(SysParams__Mx,SysParams__Mz);
   Psi3=zeros(SysParams__Mx,SysParams__Mz);
   n11=zeros(SysParams__Mx,SysParams__Mz);
   n12=zeros(SysParams__Mx,SysParams__Mz);
   n13=zeros(SysParams__Mx,SysParams__Mz);
   n22=zeros(SysParams__Mx,SysParams__Mz);
   n23=zeros(SysParams__Mx,SysParams__Mz);
   n33=zeros(SysParams__Mx,SysParams__Mz);
   nTOF11=zeros(SysParams__Mx,SysParams__Mz);
   nTOF12=zeros(SysParams__Mx,SysParams__Mz);
   nTOF13=zeros(SysParams__Mx,SysParams__Mz);
   nTOF22=zeros(SysParams__Mx,SysParams__Mz);
   nTOF23=zeros(SysParams__Mx,SysParams__Mz);
   nTOF33=zeros(SysParams__Mx,SysParams__Mz);
   tmp_Mx=zeros(SysParams__Mx,SysParams__Mz);
   
   % Fill arrays
   for j=1:24
       for jx=1:SysParams__Mx
          fseek(fileID,8*SysParams__Mz*((my_value-1)+SysParams__My*(jx-1+(j-1)*SysParams__Mx)),'bof');
          tmp_array=fread(fileID,SysParams__Mz,'double=>double'); 
          tmp_Mx(jx,:)=tmp_array.';
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
   clear tmp_array tmp_Mx;
   
   % Filling up x and z
   x=((1:SysParams__Mx)-0.5*(1+SysParams__Mx))*SysParams__ax;
   z=((1:SysParams__Mz)-0.5*(1+SysParams__Mz))*SysParams__az;
   
   % Hedgehogize
   if (Hedgehogize==1 || Hedgehogize==2)
       X=ones(SysParams__Mx,SysParams__Mz);
       Y=ones(SysParams__Mx,SysParams__Mz);
       Z=ones(SysParams__Mx,SysParams__Mz);
       for i1=1:SysParams__Mx
           X(i1,:)=X(i1,:)*x(i1);
       end
       for i1=1:SysParams__Mz
           Z(:,i1)=Z(:,i1)*z(i1);
       end
       Y=Y*y_value;
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