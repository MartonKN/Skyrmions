function plotLine1(lineVector,ArrayToPlot,figureNum,figureTitle)
     % Deleting previous version of the figure;
     fg=figure(figureNum);
     delete(fg);
     figure(figureNum);
     
     plot(lineVector,ArrayToPlot,'r-');
     title(figureTitle);
end