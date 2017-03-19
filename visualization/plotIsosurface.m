function plotIsosurface(xPlot,yPlot,zPlot,ArrayToPlot,isoValue, ...
                       figureNum,figureTitle,xLabel,yLabel,zLabel)
     % ArrayToPlot has to be reshaped because of the idiotic data structure
     % of the MATLAB plot functions. Therefore the 1st and 2nd coordinate
     % of ArrayToPlot has to be interchanged.
     s=size(ArrayToPlot);
     tmp=zeros(s(2),s(1),s(3));
     for j1=1:s(1)
         for j2=1:s(2)
             for j3=1:s(3)
                 tmp(j2,j1,j3)=ArrayToPlot(j1,j2,j3);
             end
         end
     end
     ArrayToPlot=tmp;
     clear tmp;
     
     % Deleting previous version of the figure;
     fg=figure(figureNum);
     delete(fg);
     figure(figureNum);
     
     % Plotting ArrayToPlot
     sPlotArea=subplot(1,1,1); 
     set(sPlotArea,'Units','normalized');
     set(sPlotArea,'Position',[0.1,0.1,0.8,0.8]);
     [XPlot,YPlot,ZPlot]=meshgrid(xPlot,yPlot,zPlot);
     p=patch(isosurface(XPlot,YPlot,ZPlot,ArrayToPlot,isoValue)); 
     isonormals(XPlot,YPlot,ZPlot,ArrayToPlot,p); 
     set(p,'FaceColor','red','EdgeColor','none'); 
     daspect([1,1,1]); 
     view(3); 
     camlight; 
     lighting phong; 
     grid on;
     
     title(figureTitle);
     xlabel(xLabel);
     ylabel(yLabel);
     zlabel(zLabel);
end