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
[TT,gentaExtra,gentaIntra,azithroExtra, azithroIntra, timeRunVal,combinationMIC, Loewe_ExtraVal,Loewe_IntraValC, potent_ex_val, potent_in_val]=evalc('dualTreatmentGentaAndAzithro(LHSsamplesTreatment,clearedOutcomes_untreatedModel);');
toc