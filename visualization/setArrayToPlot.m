function ArrayToPlot=setArrayToPlot(Psi1,Psi2,Psi3,n11,n12,n13,n22,n23,n33,...
                                    nTOF11,nTOF12,nTOF13,nTOF22,nTOF23,nTOF33,...
                                    WhatToPlot,Component,whichBasis)
    % According what the user sets in the GUI, different arrays have
    % to be plotted:
    % 'WhatToPlot'   -->     'ArrayToPlot'
    % Abs(Psi)               abs(Psi)
    % Arg(Psi)               angle(Psi)
    % Re(Psi)                real(Psi)
    % Im(Psi)                imag(Psi)
    % n                      n
    % nTOF                   nTOF
    % Abs(Psi)^2             abs(Psi)^2
    %
    % 'Component'
    % 1
    % 2
    % 3
    % 1+2+3
    % |(1,2,3)|
    % 1-2
    % 2-3
    % 3-1
    
    % whichBasis:
    %    ==1: hedgehog basis (used in the GPE solver C++ code)
    %    ==2: eigenbasis of the z spinoperator
    % In order to transform the quantities from the hedgehog basis to the
    % spin basis, we introduce the matrix 
    trfMx=[1,    0,       1 ;...
           1i,   0,      -1i;...
           0,    sqrt(2), 0 ]'/sqrt(2);
    % With this we transform the quantities as:
    %  <<hedgehog basis>>           <<spin basis>>
    %         Psi                   trfMx*Psi
    %          n                    trfMx'.'*n*trfMx.'
    %         nTOF                  trfMx'.'*nTOF*trfMx.'
    
    % NOTE: This function also deletes Psi-s in order to avoid
    % superfluous memory consumption.
    
    
    
    i=sqrt(-1);
    
    if strcmpi(WhatToPlot,'abs(Psi)') || ...
       strcmpi(WhatToPlot,'arg(Psi)') || ...
       strcmpi(WhatToPlot,'re(Psi)')  || ...
       strcmpi(WhatToPlot,'im(Psi)')  || ...
       strcmpi(WhatToPlot,'abs(Psi)^2')
        clear n11 n12 n13 n22 n23 n33 nTOF11 nTOF12 nTOF13 nTOF22 nTOF23 nTOF33;
        
        if whichBasis==2
            plotArray1=trfMx(1,1)*Psi1+trfMx(1,2)*Psi2+trfMx(1,3)*Psi3;
            plotArray2=trfMx(2,1)*Psi1+trfMx(2,2)*Psi2+trfMx(2,3)*Psi3;
            plotArray3=trfMx(3,1)*Psi1+trfMx(3,2)*Psi2+trfMx(3,3)*Psi3;
        else
            plotArray1=Psi1;
            plotArray2=Psi2;
            plotArray3=Psi3;
        end
        clear Psi1 Psi2 Psi3;
    elseif strcmpi(WhatToPlot,'n')
        clear Psi1 Psi2 Psi3 nTOF11 nTOF12 nTOF13 nTOF22 nTOF23 nTOF33;

        if whichBasis==2
            % Transformation corresponding to trfMx.
            plotArray1=zeros(size(n11));
            plotArray2=zeros(size(n11));
            plotArray3=zeros(size(n11));
            n21=conj(n12); n31=conj(n13); n32=conj(n23);
            for j1=1:3
                for j2=1:3
                    for j3=1:3
                        eval(['plotArray',num2str(j1),'=plotArray',num2str(j1),...
                              '+conj(trfMx(',num2str(j1),',',num2str(j2),'))*n',...
                              num2str(j2),num2str(j3),'*trfMx(',num2str(j1),',',...
                              num2str(j3),');']);
                    end
                end
            end
        else
            plotArray1=n11;
            plotArray2=n22;
            plotArray3=n33;
        end
        clear n11 n12 n13 n22 n23 n33;
    elseif strcmpi(WhatToPlot,'ntof')
        clear Psi1 Psi2 Psi3 n11 n12 n13 n22 n23 n33;

        if whichBasis==2
            % Transformation corresponding to trfMx.
            plotArray1=zeros(size(nTOF11));
            plotArray2=zeros(size(nTOF11));
            plotArray3=zeros(size(nTOF11));
            nTOF21=conj(nTOF12); nTOF31=conj(nTOF13); nTOF32=conj(nTOF23);
            for j1=1:3
                for j2=1:3
                    for j3=1:3
                        eval(['plotArray',num2str(j1),'=plotArray',num2str(j1),...
                              '+conj(trfMx(',num2str(j1),',',num2str(j2),'))*nTOF',...
                              num2str(j2),num2str(j3),'*trfMx(',num2str(j1),',',...
                              num2str(j3),');']);
                    end
                end
            end
        else
            plotArray1=nTOF11;
            plotArray2=nTOF22;
            plotArray3=nTOF33;
        end
        clear nTOF11 nTOF12 nTOF13 nTOF22 nTOF23 nTOF33;
    else
        clear Psi1 Psi2 Psi3 n11 n12 n13 n22 n23 n33 ...
              nTOF11 nTOF12 nTOF13 nTOF22 nTOF23 nTOF33;
    end
    
    switch lower(WhatToPlot)
        case 'abs(psi)'
            plotArray1=abs(plotArray1);
            plotArray2=abs(plotArray2);
            plotArray3=abs(plotArray3);
        case 'arg(psi)'
            plotArray1=angle(plotArray1);
            plotArray2=angle(plotArray2);
            plotArray3=angle(plotArray3);
        case 're(psi)'
            plotArray1=real(plotArray1);
            plotArray2=real(plotArray2);
            plotArray3=real(plotArray3);
        case 'im(psi)'
            plotArray1=imag(plotArray1);
            plotArray2=imag(plotArray2);
            plotArray3=imag(plotArray3);
        case 'n'
        case 'ntof'
        case 'abs(psi)^2'
            plotArray1=abs(plotArray1).^2;
            plotArray2=abs(plotArray2).^2;
            plotArray3=abs(plotArray3).^2;
        otherwise
            ArrayToPlot=[];
            disp(['Inappropriate value of parameter WhatToPlot in function setArrayToPlot: ', WhatToPlot]);
    end
    
    Component=lower(strtrim(Component));
    if strcmp(Component,'1')
        ArrayToPlot=plotArray1;
    elseif strcmp(Component,'2')
        ArrayToPlot=plotArray2;
    elseif strcmp(Component,'3')
        ArrayToPlot=plotArray3;
    elseif strcmp(Component,'1+2+3')
        ArrayToPlot=plotArray1+plotArray2+plotArray3;
    elseif strcmp(Component,'|(1,2,3)|')
        ArrayToPlot=sqrt(abs(plotArray1).^2+abs(plotArray2).^2+abs(plotArray3).^2);
    elseif strcmp(Component,'1-2')
        ArrayToPlot=plotArray1-plotArray2;
    elseif strcmp(Component,'2-3')
        ArrayToPlot=plotArray2-plotArray3;
    elseif strcmp(Component,'3-1')
        ArrayToPlot=plotArray3-plotArray1;
    else
        ArrayToPlot=[];
        disp(['Inappropriate value of parameter Component in function setArrayToPlot: ', Component]);
    end
    clear plotArray1 plotArray2 plotArray3;
end