function [CmaxAboveMIC,summary_Cmax_above_MIC]=Cmax(antibiotic,~)
%%Cmax/MIC
numSamples=size(antibiotic,2);

MICchecked=1;
CmaxAboveMIC=NaN(length(MICchecked),numSamples);
summary_Cmax_above_MIC=NaN(length(MICchecked),3);

for i=1:length(MICchecked)
    for j=1:numSamples
        
        CmaxAboveMIC(i,j)=max(antibiotic(:,j))/MICchecked(i);
        
    end
end

for i=1:length(MICchecked)
    summary_Cmax_above_MIC(i,:)=quantile(CmaxAboveMIC(i,:),[0.5,0.025,0.975]);
end
end