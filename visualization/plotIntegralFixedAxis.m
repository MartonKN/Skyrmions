function plotIntegralFixedAxis(whichAxis,DataArray,plotDim,figureNum,colorMap,figureTitle,filename_SysParams)
load(filename_SysParams);

x=((-0.5*(SysParams__Mx-1)):1:(0.5*(SysParams__Mx-1)))*SysParams__ax;
y=((-0.5*(SysParams__My-1)):1:(0.5*(SysParams__My-1)))*SysParams__ay;
z=((-0.5*(SysParams__Mz-1)):1:(0.5*(SysParams__Mz-1)))*SysParams__az;

if whichAxis==1
    plotArray=zeros(SysParams__My,SysParams__Mz);
    for jy=1:SysParams__My
        for jz=1:SysParams__Mz
            plotArray(jy,jz)=sum(DataArray(:,jy,jz));
        end
    end
    label1='y'; label2='z';
    plotVec1=y; plotVec2=z;
elseif whichAxis==2
    plotArray=zeros(SysParams__Mx,SysParams__Mz);
    for jx=1:SysParams__Mx
        for jz=1:SysParams__Mz
            plotArray(jx,jz)=sum(DataArray(jx,:,jz));
        end
    end
    label1='x'; label2='z';
    plotVec1=x; plotVec2=z;
elseif whichAxis==3
    plotArray=zeros(SysParams__Mx,SysParams__My);
    for jx=1:SysParams__Mx
        for jy=1:SysParams__My
            plotArray(jx,jy)=sum(DataArray(jx,jy,:));
        end
    end
    label1='x'; label2='y';
    plotVec1=x; plotVec2=y;
end


plotFlat(plotVec1,plotVec2,plotArray,plotDim,figureNum,colorMap,...
         figureTitle,label1,label2);
end