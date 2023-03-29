function [areaAboveMIC,summary_area_above_MIC]=tot_area_under_curve_multiple_doses(antibiotic,~, numDoses,timeBetweenDose,LHS_treatment)
%calculating the total AUC 0 to infinity and dividing by the MIC

numSamples=size(antibiotic,2);
MICchecked=1;
totareaAboveMIC=NaN(1,numSamples);
areaAboveMIC=NaN(length(MICchecked),numSamples);
summary_area_above_MIC=NaN(length(MICchecked),3);

for j=1:numSamples
    cumulativeArea=0;
    
    a=find(antibiotic(:,j),1);%in the cocentration curve the first index point that drug is adminsitered.
    b=a+timeBetweenDose;%the last index before the next dose is administered
    
    for k=1:numDoses
        
        if k==numDoses%last doing interval so need to integrate upto infinity
            
            C0=antibiotic(a,j);
            delta=LHS_treatment(j,6);
            fun = @(x)C0*exp(-delta*x);
            t=0;
            totalAreaBelow= integral(fun,t,Inf);
            cumulativeArea= cumulativeArea+totalAreaBelow;
        else
            
            totalAreaBelow=trapz(antibiotic(a:b,j));%gives the total area under the curve
            cumulativeArea= cumulativeArea+totalAreaBelow;
            a=b;
            b=a+timeBetweenDose;
        end
        
    end
    
    totareaAboveMIC(1,j)=cumulativeArea;
  
end


for j=1:numSamples
    for i=1:length(MICchecked)
        areaAboveMIC(i,j)=totareaAboveMIC(1,j)./MICchecked(i);
    end
    
end

for i=1:length(MICchecked)
    summary_area_above_MIC(i,:)=quantile(areaAboveMIC(i,:),[0.5,0.025,0.975]);
end

end