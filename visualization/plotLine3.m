function plotLine3(lineVector,ArrayToPlot1,legend1,...
                   ArrayToPlot2,legend2,ArrayToPlot3,legend3,...
                   figureNum,figureTitle)
     % Deleting previous version of the figure;
     fg=figure(figureNum);
     delete(fg);
     figure(figureNum);
     
     plot(lineVector,ArrayToPlot1,'r-',lineVector,ArrayToPlot2,'g-',...
          lineVector,ArrayToPlot3,'b-');
     legend(legend1,legend2,legend3);
     title(figureTitle);
end