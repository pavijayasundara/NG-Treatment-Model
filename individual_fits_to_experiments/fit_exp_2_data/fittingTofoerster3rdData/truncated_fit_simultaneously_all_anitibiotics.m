clear
close all
numAntibioitcs=9;
MICest=NaN(numAntibioitcs,1);
Hill_k=NaN(numAntibioitcs,1);
errors=NaN(numAntibioitcs,1);
pi_min=NaN(numAntibioitcs,1);
confIntMIC=NaN(numAntibioitcs,2);
confIntHill_k=NaN(numAntibioitcs,2);
confIntphi_min=NaN(numAntibioitcs,2);
phi_max=phi_max_all_antibioitcs();
concVals=NaN(12,numAntibioitcs);
x1=2;
x2=3;
figure('units','normalized','outerposition',[0 0 1 1])

subplot(x1,x2,1)
load('C:\Users\z5137734\OneDrive - UNSW\PhD\treatment model\treatment model with new equations\Foerester model reproduce\fittingMyModelToFoerster\data_3\gentamicin\gentamicin_data')
[fittedVals,errorVals,confint,SSEconc]= truncated_fit_simultaneouslyGenta(NGVals,antibioticConcVals,'gentamicin',phi_max(1),-15);
MICest(1)=fittedVals(1);
Hill_k(1)=fittedVals(2);
pi_min(1)=fittedVals(3);
errors(1)=errorVals;
confIntMIC(1,1:2)=confint(1,1:2);
confIntHill_k(1,1:2)=confint(2,1:2);
confIntphi_min(1,1:2)=confint(3,1:2);
concVals(:,1)=antibioticConcVals;
% SSE_conc1=NaN(8,2);%errors at each concnetration. 11 concentrations 
% SSE_conc1(:,1)=antibioticConcVals(2:end-3);
% SSE_conc1(:,2)=SSEconc;


subplot(x1,x2,2)
%clear
load('C:\Users\z5137734\OneDrive - UNSW\PhD\treatment model\treatment model with new equations\Foerester model reproduce\fittingMyModelToFoerster\data_3\spec\spec_data')
[fittedVals,errorVals,confint,SSEconc]=truncated_fit_simultaneouslyGenta(NGVals,antibioticConcVals,'spectinomycin',phi_max(2),-8);
MICest(2)=fittedVals(1);
Hill_k(2)=fittedVals(2);
errors(2)=errorVals;
pi_min(2)=fittedVals(3);
confIntMIC(2,1:2)=confint(1,1:2);
confIntHill_k(2,1:2)=confint(2,1:2);
confIntphi_min(2,1:2)=confint(3,1:2);
concVals(:,2)=antibioticConcVals;

subplot(x1,x2,3)
%clear
load('C:\Users\z5137734\OneDrive - UNSW\PhD\treatment model\treatment model with new equations\Foerester model reproduce\fittingMyModelToFoerster\data_3\azithromycin\azithroData')
[fittedVals,errorVals,confint,SSEconc]= truncated_fit_simultaneously(NGVals,antibioticConcVals,'azithromycin',phi_max(3),-3.3);
MICest(3)=fittedVals(1);
Hill_k(3)=fittedVals(2);
errors(3)=errorVals;
pi_min(3)=fittedVals(3);
confIntMIC(3,1:2)=confint(1,1:2);
confIntHill_k(3,1:2)=confint(2,1:2);
confIntphi_min(3,1:2)=confint(3,1:2);
% SSE_conc3=NaN(11,2);%errors at each concnetration. 11 concentrations 
% SSE_conc3(:,1)=antibioticConcVals(2:end);
% SSE_conc3(:,2)=SSEconc;
concVals(:,3)=antibioticConcVals;

subplot(x1,x2,4)
%clear
load('C:\Users\z5137734\OneDrive - UNSW\PhD\treatment model\treatment model with new equations\Foerester model reproduce\fittingMyModelToFoerster\data_3\penicillin\penicillin_data')
[fittedVals,errorVals,confint,SSEconc]= truncated_fit_simultaneouslyPenicillin(NGVals,antibioticConcVals,'penicillin',phi_max(4),-3.6);
MICest(4)=fittedVals(1);
Hill_k(4)=fittedVals(2);
errors(4)=errorVals;
pi_min(4)=fittedVals(3);
confIntMIC(4,1:2)=confint(1,1:2);
confIntHill_k(4,1:2)=confint(2,1:2);
confIntphi_min(4,1:2)=confint(3,1:2);
% SSE_conc4=NaN(11,2);%errors at each concnetration. 11 concentrations 
% SSE_conc4(:,1)=antibioticConcVals(2:end);
% SSE_conc4(:,2)=SSEconc;
concVals(1:10,4)=antibioticConcVals;

subplot(x1,x2,5)
%clear
load('C:\Users\z5137734\OneDrive - UNSW\PhD\treatment model\treatment model with new equations\Foerester model reproduce\fittingMyModelToFoerster\data_3\ceftriaxone\ceftriaxone_data')
[fittedVals,errorVals,confint,SSEconc]= truncated_fit_simultaneouslyBetaLactams(NGVals,antibioticConcVals,'ceftriaxone',phi_max(5),-1.2);
MICest(5)=fittedVals(1);
Hill_k(5)=fittedVals(2);
errors(5)=errorVals;
pi_min(5)=fittedVals(3);
confIntMIC(5,1:2)=confint(1,1:2);
confIntHill_k(5,1:2)=confint(2,1:2);
confIntphi_min(5,1:2)=confint(3,1:2);
% SSE_conc5=NaN(11,2);%errors at each concnetration. 11 concentrations 
% SSE_conc5(:,1)=antibioticConcVals(2:end);
% SSE_conc5(:,2)=SSEconc;
concVals(1:9,5)=antibioticConcVals;


subplot(x1,x2,6)
%clear
load('C:\Users\z5137734\OneDrive - UNSW\PhD\treatment model\treatment model with new equations\Foerester model reproduce\fittingMyModelToFoerster\data_3\cefixime\cefixime_data')
[fittedVals,errorVals,confint,SSEconc]= truncated_fit_simultaneouslyBetaLactams(NGVals,antibioticConcVals,'cefixime',phi_max(6),-1.5);
MICest(6)=fittedVals(1);
Hill_k(6)=fittedVals(2);
errors(6)=errorVals;
pi_min(6)=fittedVals(3);
confIntMIC(6,1:2)=confint(1,1:2);
confIntHill_k(6,1:2)=confint(2,1:2);
confIntphi_min(6,1:2)=confint(3,1:2);
concVals(1:9,6)=antibioticConcVals;

x1=1;
x2=3;
figure('units','normalized','outerposition',[0 0 1 1])
subplot(x1,x2,1)
%clear
load('C:\Users\z5137734\OneDrive - UNSW\PhD\treatment model\treatment model with new equations\Foerester model reproduce\fittingMyModelToFoerster\data_3\chloro\chloramphenicol_data')
[fittedVals,errorVals,confint,SSEconc]= truncated_fit_simultaneouslyChloramphenicol(NGVals,antibioticConcVals,'chloramphenicol',phi_max(7),-0.5);
MICest(7)=fittedVals(1);
Hill_k(7)=fittedVals(2);
errors(7)=errorVals;
pi_min(7)=fittedVals(3);
confIntMIC(7,1:2)=confint(1,1:2);
confIntHill_k(7,1:2)=confint(2,1:2);
confIntphi_min(7,1:2)=confint(3,1:2);
concVals(1:11,7)=antibioticConcVals;


subplot(x1,x2,2)
%clear
load('C:\Users\z5137734\OneDrive - UNSW\PhD\treatment model\treatment model with new equations\Foerester model reproduce\fittingMyModelToFoerster\data_3\tetra\tetracycline_data')
[fittedVals,errorVals,confint,SSEconc]= truncated_fit_simultaneously(NGVals,antibioticConcVals,'tetracycline',phi_max(8),-0.0005);
MICest(8)=fittedVals(1);
Hill_k(8)=fittedVals(2);
errors(8)=errorVals;
pi_min(8)=fittedVals(3);
confIntMIC(8,1:2)=confint(1,1:2);
confIntHill_k(8,1:2)=confint(2,1:2);
confIntphi_min(8,1:2)=confint(3,1:2);
concVals(:,8)=antibioticConcVals;


subplot(x1,x2,3)
%clear
load('C:\Users\z5137734\OneDrive - UNSW\PhD\treatment model\treatment model with new equations\Foerester model reproduce\fittingMyModelToFoerster\data_3\cipro\ciprofloxacin_data')
[fittedVals,errorVals,confint,SSEconc]= truncated_fit_simultaneouslyCipro(NGVals,antibioticConcVals,'ciprofloxacin',phi_max(9),-1);
MICest(9)=fittedVals(1);
Hill_k(9)=fittedVals(2);
errors(9)=errorVals;
pi_min(9)=fittedVals(3);
confIntMIC(9,1:2)=confint(1,1:2);
confIntHill_k(9,1:2)=confint(2,1:2);
confIntphi_min(9,1:2)=confint(3,1:2);
concVals(:,9)=antibioticConcVals;

results=[pi_min,confIntphi_min,MICest,confIntMIC,Hill_k,confIntHill_k,errors];
%%%errors=log10(errors);