function [timeAboveMIC,summary_t_above_MIC]=time_above_MIC_multiple_doses(antibiotic,timeRun, numDoses,timeBetweenDose)

numSamples=size(antibiotic,2);
MICchecked=1;
timeAboveMIC=NaN(length(MICchecked),numSamples);
summary_t_above_MIC=NaN(length(MICchecked),3);


for i=1:length(MICchecked)
    for j=1:numSamples
        cumulativeTIme=0;
               
        a=find(antibiotic(:,j),1);%in the cocentration curve the first index point that drug is adminsitered. 
        b=a+timeBetweenDose; %the last index before the next dose is administered
        
        for k=1:numDoses
            if max(antibiotic(a:b,j))<MICchecked(i)
                cumulativeTIme=0;
            else
            [d,ix] = min(abs(antibiotic(a:b,j)-MICchecked(i)));
            timeAboveInthisSegment=timeRun(ix+a,j)-timeRun(a);
            cumulativeTIme= cumulativeTIme+timeAboveInthisSegment;
            a=b;
            if k==numDoses-1
               b=timeRun(end); 
            else
                b=a+timeBetweenDose;
            end
            end
        end
        
        timeAboveMIC(i,j)=cumulativeTIme;
        
    end
end

for i=1:length(MICchecked)
    summary_t_above_MIC(i,:)=quantile(timeAboveMIC(i,:),[0.5,0.025,0.975]);
end
end