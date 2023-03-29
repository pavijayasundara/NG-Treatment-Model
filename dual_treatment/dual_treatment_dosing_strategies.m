clear
close all
tic

load('treatementParaSamples')
LHS_output=LHS_output(2:end,:);
LHS_treatment=cell2mat(LHS_output);
load('LHS_samples_genta')
LHS_treatment(:,[3,5,4])=LHS_genta(:,2:end);
load('LHS_samples_azithro')
LHS_treatment(:,[11,13,12])=LHS_azithro(:,2:end);
load('cleared samples untreated infection')
LHSsamplesTreatment(:,15:32)=LHS_treatment;

[TT,clearence_Times,gentaExtra,timeRunVal,total_BacteriaVal,combinationMIC]=evalc('dualTreatmentGentaAndAzithro(LHSsamplesTreatment,clearedOutcomes_untreatedModel);');
toc