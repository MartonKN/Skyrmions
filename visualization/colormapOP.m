function cMap=colormapOP()
cMap=jet(500); 
cMap=cMap(floor(length(cMap)*0.3):length(cMap),:); 
x=linspace(0,1,floor(length(cMap)/3))';
cMap=[cMap;[0.3*x+cMap(length(cMap),1)*(1-x),cMap(length(cMap),2)*(1-x),cMap(length(cMap),3)*(1-x)]];
end