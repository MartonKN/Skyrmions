function plotEllipsoid(XPlot,YPlot,ZPlot,ArrayToPlot,figureNum,colorMap,...
                       figureTitle,xLabel,yLabel,zLabel)
     % Deleting previous version of the figure;
     fg=figure(figureNum);
     delete(fg);
     figure(figureNum);
     
     % Plotting color code
     minValue=min(ArrayToPlot); minValue=min(minValue);
     maxValue=max(ArrayToPlot); maxValue=max(maxValue);
     if (minValue==maxValue)
         minValue=minValue-1e-5;
         maxValue=maxValue+1e-5;
     end
     sColorCode=subplot(1,2,2);
     set(sColorCode,'Units','normalized');
     set(sColorCode,'Position',[0.92,0.115,0.06,0.82]);
     xColorCode=[0,1];
     yColorCode=linspace(minValue,maxValue,256);
     [XColorCode,YColorCode]=meshgrid(xColorCode,yColorCode); 
     if(strcmpi(colorMap,'vortices'))
         [colorMapPlotting,colorMapColorCode]=...
             setVorticesColormap(minValue,maxValue);
         colormap(colorMapColorCode);
     else
         colormap(colorMap); 
     end
     pcolor(XColorCode,YColorCode,YColorCode); 
     view(0,90); 
     shading 'flat';
     
     % Plotting ArrayToPlot
     sPlotArea=subplot(1,2,1); 
     set(sPlotArea,'Units','normalized');
     set(sPlotArea,'Position',[0.12,0.15,0.75,0.75]);
     if(strcmpi(colorMap,'vortices'))
         colormap(colorMapPlotting);
     else
         colormap(colorMap); 
     end
     surf(XPlot,YPlot,ZPlot,ArrayToPlot);
     title(figureTitle);
     xlabel(xLabel);
     ylabel(yLabel);
     zlabel(zLabel);
     shading 'flat';
end