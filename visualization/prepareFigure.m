function prepareFigure(figureNum,colorMap,minValue,maxValue)
     % Prepares the figure for plotting:
     % Calls figure(figureNum); and separates one area for the color scale
     % on the right. minValue and maxValue stand for the minimal and 
     % maximal value of the function to be plotted. 
     
     fg=figure(figureNum);
     delete(fg);
     figure(figureNum);
     sColorCode=subplot(1,2,2);
     set(sColorCode,'Units','normalized');
     set(sColorCode,'Position',[0.92,0.1,0.06,0.8]);
     xColorCode=[0,1]; 
     if (minValue==maxValue)
         minValue=minValue-1e-5;
         maxValue=maxValue+1e-5;
     end
     yColorCode=linspace(minValue,maxValue,256); 
     [X,Y]=meshgrid(xColorCode,yColorCode); 
     colormap(colorMap); 
     pcolor(X,Y,Y); 
     view(0,90); 
     shading 'flat';
     sPlotArea=subplot(1,2,1); 
     set(sPlotArea,'Units','normalized');
     set(sPlotArea,'Position',[0.09,0.1,0.75,0.8]);     
end