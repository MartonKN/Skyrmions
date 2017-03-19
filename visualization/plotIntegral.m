function plotIntegral(theta,phi,DataArray,plotDim,figureNum,colorMap,figureTitle,filename_SysParams)
load(filename_SysParams);

transfMx=[cos(phi),sin(phi),0;-sin(phi),cos(phi),0;0,0,1]*...
         [1,0,0;0,cos(theta),sin(theta);0,-sin(theta),cos(theta)];
x=((-0.5*(SysParams__Mx-1)):1:(0.5*(SysParams__Mx-1)))*SysParams__ax;
y=((-0.5*(SysParams__My-1)):1:(0.5*(SysParams__My-1)))*SysParams__ay;
z=((-0.5*(SysParams__Mz-1)):1:(0.5*(SysParams__Mz-1)))*SysParams__az;

xMin=-0.5*(SysParams__Mx-1)*SysParams__ax;
xMax= 0.5*(SysParams__Mx-1)*SysParams__ax;
yMin=-0.5*(SysParams__My-1)*SysParams__ay;
yMax= 0.5*(SysParams__My-1)*SysParams__ay;
zMin=-0.5*(SysParams__Mz-1)*SysParams__az;
zMax= 0.5*(SysParams__Mz-1)*SysParams__az;

edges(:,1)=[xMin;yMin;zMin];
edges(:,2)=[xMin;yMin;zMax];
edges(:,3)=[xMin;yMax;zMin];
edges(:,4)=[xMin;yMax;zMax];
edges(:,5)=[xMax;yMin;zMin];
edges(:,6)=[xMax;yMin;zMax];
edges(:,7)=[xMax;yMax;zMin];
edges(:,8)=[xMax;yMax;zMax];

edges=transfMx*edges;
newXMin=min(edges(1,:));
newXMax=max(edges(1,:));
newYMin=min(edges(2,:));
newYMax=max(edges(2,:));

newXArray=linspace(newXMin,newXMax,SysParams__Mx);
newYArray=linspace(newYMin,newYMax,SysParams__My);
plotArray=zeros(SysParams__Mx,SysParams__My);
newX_stepsize=(newXMax-newXMin)/(SysParams__Mx-1);
newY_stepsize=(newYMax-newYMin)/(SysParams__My-1);

% Dilute data
sX=max([1,round(SysParams__Mx/32)]);
sY=max([1,round(SysParams__My/32)]);
sZ=max([1,round(SysParams__Mz/32)]);

for jx=1:sX:SysParams__Mx
    for jy=1:sY:SysParams__My
        for jz=1:sZ:SysParams__Mz
            r=[x(jx);y(jy);z(jz)];
            r=transfMx*r;
            index1=floor((r(1)-newXMin)/newX_stepsize)+1;
            index2=floor((r(2)-newYMin)/newY_stepsize)+1;
            plotArray(index1,index2)=plotArray(index1,index2)+DataArray(jx,jy,jz);
        end
    end
end

plotFlat(newXArray,newYArray,plotArray,plotDim,figureNum,colorMap,...
                  figureTitle,'','');
end