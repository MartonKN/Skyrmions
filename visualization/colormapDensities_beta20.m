function cMap=colormapDensities_beta20()
% Jet / green->white
cMap=jet(500); %cMap=0.9*cMap+0.1*ones(size(cMap));
x=linspace(0,1,floor(length(cMap)/20))';

%   Add dark colors to the ends of the jet colormap
tmp=[x*cMap(1,1),x*cMap(1,2),0.1*(1-x)+x*cMap(1,3)];
cMap=[tmp;cMap];
tmp=[0.1*x+(1-x)*cMap(length(cMap),1),(1-x)*cMap(length(cMap),2),(1-x)*cMap(length(cMap),3)];
cMap=[cMap;tmp];

%   Change the middle of the colorbar to white
tmp=cMap(floor(length(cMap)*0.4):floor(length(cMap)*0.5),:);
x=linspace(0,1,length(tmp))';
tmp=[x+(1-x)*tmp(1,1),x+(1-x)*tmp(1,2),x+(1-x)*tmp(1,3)];
cMap(floor(length(cMap)*0.4):floor(length(cMap)*0.5),:)=tmp;

tmp=cMap(floor(length(cMap)*0.5):floor(length(cMap)*0.6),:);
x=linspace(0,1,length(tmp))';
tmp=[1-x+(x)*tmp(length(tmp),1),1-x+(x)*tmp(length(tmp),2),1-x+(x)*tmp(length(tmp),3)];
cMap(floor(length(cMap)*0.5):floor(length(cMap)*0.6),:)=tmp;



% % Jet + white + black
% % negative part
% cMap1=jet(500); 
% cMap1=cMap1(1:floor(length(cMap1)*0.4),:); 
% x=linspace(0,1,floor(length(cMap1)/2.5))';  
% y=linspace(0,1,floor(length(cMap1)/2.5))'; 
% tmp=[x*cMap1(1,1),x*cMap1(1,2),0.2*(1-x)+x*cMap1(1,3)];
% cMap1=[tmp;cMap1]; 
% cMap1=flipud(cMap1); 
% tmp=[1-y+y*cMap1(1,1),1-y+y*cMap1(1,2),1-y+y*cMap1(1,3)]; 
% cMap1=[tmp;cMap1];
% % positive part
% cMap2=jet(500); 
% cMap2=cMap2(floor(length(cMap2)*0.6:length(cMap2)),:); 
% cMap2=flipud(cMap2);
% x=linspace(0,1,floor(length(cMap2)/3))';  
% y=linspace(0,1,floor(length(cMap2)/2))'; 
% tmp=[x*cMap2(1,1),x*cMap2(1,2),x*cMap2(1,3)];
% cMap2=[tmp;cMap2]; 
% cMap2=flipud(cMap2); 
% tmp=[1-y+y*cMap2(1,1),1-y+y*cMap2(1,2),1-y+y*cMap2(1,3)]; 
% cMap2=[tmp;cMap2];
% cMap2=flipud(cMap2);
% % putting the two together
% cMap=flipud([cMap2;cMap1]);


% % Jet positve values + white
% cMap=jet(500);
% cMap=cMap(floor(0.58*length(cMap)):length(cMap),:);
% x=linspace(0,1,length(cMap)/5)';
% tmp=[1-x+x*cMap(1,1),1-x+x*cMap(1,2),1-x+x*cMap(1,3)];
% cMap=[tmp;cMap];


% % Yellow, red and balck - from Pascu
% white=[1,1,1];
% black=[0,0,0];
% blue=[0,0,1];
% green=[0,1,0];
% yellow=[1,1,0]; 
% red=[1,0,0]; 
% 
% c(1,:)=white;
% c(2,:)=yellow;
% c(3,:)=red;
% c(4,:)=black;
% 
% colorLength=[40,70,70,100];
% x=linspace(0,1,max(colorLength))';
% cMap=[];
% 
% for j1=1:3
%     tmpCMap=[(1-linspace(0,1,colorLength(j1))')*c(j1,1)+...
%              linspace(0,1,colorLength(j1))'*c(j1+1,1),...
%              (1-linspace(0,1,colorLength(j1))')*c(j1,2)+...
%              linspace(0,1,colorLength(j1))'*c(j1+1,2),...
%              (1-linspace(0,1,colorLength(j1))')*c(j1,3)+...
%              linspace(0,1,colorLength(j1))'*c(j1+1,3)];
%     cMap=[cMap;tmpCMap];
% end

% %Greiner's color scale
% white=[1,1,1]; 
% blue=[0,0,1]; 
% green=[0,1,0];
% yellow=[1,1,0]; 
% red=[1,0,0]; 
% 
% c(1,:)=[0,0,0.3];
% c(2,:)=[0,0,0.6];
% c(3,:)=blue;
% c(4,:)=green;
% c(5,:)=white;
% c(6,:)=yellow;
% c(7,:)=red;
% c(8,:)=[0.6,0,0];
% c(9,:)=[0.3,0,0];
% 
% x=linspace(0,1,100)'; 
% cMap=[];
% 
% for j1=1:8
%     tmpCMap=[(1-x)*c(j1,1)+x*c(j1+1,1),...
%                  (1-x)*c(j1,2)+x*c(j1+1,2),...
%                  (1-x)*c(j1,3)+x*c(j1+1,3)];
%     cMap=[cMap;tmpCMap];
% end

end