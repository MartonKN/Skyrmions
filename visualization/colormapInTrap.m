function cMap=colormapInTrap()
white=[1,1,1]; 
darkblue=[50,89,164]/256;
green=[0,172,88]/256;
yellow=[254,230,60]/256; 
red=[203,63,62]/256; 

c(1,:)=darkblue/1.2;
c(2,:)=white;
c(3,:)=red/1.2;

x=linspace(0,1,100)'.^1;
cMap=[];

for j1=1:2
    tmpCMap=[(1-x)*c(j1,1)+x*c(j1+1,1),...
             (1-x)*c(j1,2)+x*c(j1+1,2),...
             (1-x)*c(j1,3)+x*c(j1+1,3)];
    cMap=[cMap;tmpCMap];
end

end