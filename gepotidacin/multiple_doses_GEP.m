function [clearence_Times,timeRunVal,total_BacteriaVal,antibioticVal]=multiple_doses_GEP(LHSSamples,clearedoutcomes)
numLHSsamples=size(LHSSamples,1);

MIC_GEP=[0.5,0.125,0.25,0.5,1];
numMIC=length(MIC_GEP);
clearence_Times=NaN(numLHSsamples,numMIC);%colums are the MIC values and the rows gives the clearance times for the number of LHS samples
antibioticVal=cell(numMIC,numLHSsamples);%rows are the MIC values and the columns are the LHS
timeRunVal=cell(numMIC,numLHSsamples); %rows are the MIC values and the columns are the LHS
total_BacteriaVal=cell(numMIC,numLHSsamples); %rows are the MIC values and the columns are the LHS

initDose=1500; % initial dose that is given
for i=1:numMIC
    MICval=MIC_GEP(i);
    parfor j=1:numLHSsamples
        peakTime=clearedoutcomes(j,3); % treatmetn is intiated at the peak time of the untreated infection
        MyPara=[LHSSamples(j,:),MICval];
        [~,~,~,~,~,antibiotic,total_Bacteria,clearenceTime_with_antibiotic,timeRun]=treatment_strategy(MyPara,peakTime*24,initDose);
        if ~(isempty(clearenceTime_with_antibiotic))
            clearence_Times(j,i)=min(clearenceTime_with_antibiotic-peakTime);
        end        

        timeRunVal{i,j}=timeRun;
        total_BacteriaVal{i,j}=total_Bacteria;
        antibioticVal{i,j}=antibiotic;
    
    end
end

end
