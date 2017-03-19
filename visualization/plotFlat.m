function plotFlat(xVector,yVector,ArrayToPlot,PlotDim,figureNum,colorMap,...
                  figureTitle,xLabel,yLabel)
     % Deleting previous version of the figure;
     fg=figure(figureNum);
     delete(fg);
     figure(figureNum);
     
     % Plotting color code
     minValue=min(ArrayToPlot); minValue=min(minValue),
     maxValue=max(ArrayToPlot); maxValue=max(maxValue),
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
     colormap(colorMap); 
     pcolor(XColorCode,YColorCode,YColorCode); 
     view(0,90); 
     shading 'flat';
      
     % Plotting ArrayToPlot
     sPlotArea=subplot(1,2,1); 
     set(sPlotArea,'Units','normalized');
     set(sPlotArea,'Position',[0.11,0.115,0.70,0.82]);
     colormap(colorMap); 
     if PlotDim==2
         pcolor(xVector,yVector,ArrayToPlot');
     elseif PlotDim==3
         surf(xVector,yVector,ArrayToPlot');
     end
     title(figureTitle);
     xlabel(xLabel);
     ylabel(yLabel);
     shading 'flat';
end