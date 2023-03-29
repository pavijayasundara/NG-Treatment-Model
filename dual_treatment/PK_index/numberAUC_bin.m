function [binAUC_in_clear,binAUC_in_unclear]=numberAUC_bin(areaAboveMIC_in,clearence_Times,start,endVal,increment,finalVal,testingMIC_index)
testingMIC_index_cl=1;%%for AUC as first row shows when MIC=4
numBin=((finalVal-start)/increment)+1;
binAUC_in_clear=NaN(numBin,1);
binAUC_in_unclear=NaN(numBin,1);

a=start;
b=endVal;
for i=1:numBin
   
    if i==numBin
        indexVals=find(areaAboveMIC_in(testingMIC_index,:)>a & clearence_Times(:,testingMIC_index_cl)'<7);
    elseif i==1
        indexVals=find(areaAboveMIC_in(testingMIC_index,:)<=b & clearence_Times(:,testingMIC_index_cl)'<7); 
    else
        indexVals=find(areaAboveMIC_in(testingMIC_index,:)>a & areaAboveMIC_in(testingMIC_index,:)<=b & clearence_Times(:,testingMIC_index_cl)'<7);
    end
if any(length(indexVals))
    binAUC_in_clear(i,1)=length(indexVals);
else
    binAUC_in_clear(i,1)=NaN;
end
a=b;
b=b+increment;
end



a=start;
b=endVal;
for i=1:numBin
    
    if i==numBin
        indexVals=find(areaAboveMIC_in(testingMIC_index,:)>a & clearence_Times(:,testingMIC_index_cl)'>7);
    elseif i==1
        indexVals=find(areaAboveMIC_in(testingMIC_index,:)<=b & clearence_Times(:,testingMIC_index_cl)'>7);    
    else
        indexVals=find(areaAboveMIC_in(testingMIC_index,:)>a & areaAboveMIC_in(testingMIC_index,:)<=b & clearence_Times(:,testingMIC_index_cl)'>7);
    end

if any(length(indexVals))
    binAUC_in_unclear(i,1)=length(indexVals);
else
    binAUC_in_unclear(i,1)=NaN;
end
a=b;
b=b+increment;
end
end