function cMap=OPcolormap()

% % jet + white
% tmp=jet(100);
% c(1,:)=[1,1,1];
% c(2,:)=tmp(10,:);
% 
% x=linspace(0,1,100)'; 
% 
% cMap=[(1-x)*c(1,1)+x*c(1+1,1),...
%       (1-x)*c(1,2)+x*c(1+1,2),...
%       (1-x)*c(1,3)+x*c(1+1,3)];
% cMap=[cMap;tmp(11:100,:)];

% Red
white=[1,1,1]; 
darkblue=[50,89,164]/256; 
blue=[1,184,242]/256;
green=[0,172,88]/256;
yellow=[254,230,60]/256; 
red=[203,63,62]/256; 

c(1,:)=white;
c(2,:)=(red*1.2+yellow*0.5)/1.5;
c(3,:)=red;
c(4,:)=red/2;

x=linspace(0,1,100)'; 
cMap=[];

for j1=1:3
    tmpCMap=[(1-x)*c(j1,1)+x*c(j1+1,1),...
                 (1-x)*c(j1,2)+x*c(j1+1,2),...
                 (1-x)*c(j1,3)+x*c(j1+1,3)];
    cMap=[cMap;tmpCMap];
end

end
