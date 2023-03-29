clear
load('treatementParaSamples')% LHS samples relevant to treatment
LHS_output=LHS_output(2:end,:);
LHS_treatment=cell2mat(LHS_output);
load('cleared samples untreated infection')%natural infection parameters
LHSsamplesTreatment(:,15:20)=LHS_treatment; %append treatment LHS samples to untreated infection samples                                                                                         
[TT,clearence_Times,timeRunVal,total_BacteriaVal,...
    antibioticVal]=evalc('multiple_doses_GEP(LHSsamplesTreatment(1:5,:),clearedOutcomes_untreatedModel(1:5,:));');
