clear
close all
numAntibioitcs=9;

MICest=NaN(numAntibioitcs,1);
Hill_k=NaN(numAntibioitcs,1);
errors=NaN(numAntibioitcs,1);
pi_min=NaN(numAntibioitcs,1);

phi_max=phi_max_all_antibioitcs();

current_dir = pwd;

files = dir([current_dir '\data_exp_2\gentamicin\*.mat']);
load([files.folder '\gentamicin_data.mat'])
fittedVals=  truncated_fit_simultaneouslyGenta(NGVals,antibioticConcVals,'gentamicin',phi_max(1),-15);
MICest(1)=fittedVals(1);
Hill_k(1)=fittedVals(2);
pi_min(1)=fittedVals(3);


files = dir([current_dir '\data_exp_2\spectinomycin\*.mat']);
load([files.folder '\spec_data.mat'])
fittedVals= truncated_fit_simultaneouslyGenta(NGVals,antibioticConcVals,'spectinomycin',phi_max(2),-8);
MICest(2)=fittedVals(1);
Hill_k(2)=fittedVals(2);
pi_min(2)=fittedVals(3);

files = dir([current_dir '\data_exp_2\azithromycin\*.mat']);
load([files.folder '\azithro_data.mat'])
fittedVals= truncated_fit_simultaneously(NGVals,antibioticConcVals,'azithromycin',phi_max(3),-3.3);
MICest(3)=fittedVals(1);
Hill_k(3)=fittedVals(2);
pi_min(3)=fittedVals(3);

files = dir([current_dir '\data_exp_2\penicillin\*.mat']);
load([files.folder '\penicillin_data.mat'])
fittedVals= truncated_fit_simultaneouslyPenicillin(NGVals,antibioticConcVals,'penicillin',phi_max(4),-3.6);
MICest(4)=fittedVals(1);
Hill_k(4)=fittedVals(2);
pi_min(4)=fittedVals(3);

files = dir([current_dir '\data_exp_2\ceftriaxone\*.mat']);
load([files.folder '\ceftriaxone_data.mat'])
fittedVals= truncated_fit_simultaneouslyBetaLactams(NGVals,antibioticConcVals,'ceftriaxone',phi_max(5),-1.2);
MICest(5)=fittedVals(1);
Hill_k(5)=fittedVals(2);
pi_min(5)=fittedVals(3);



files = dir([current_dir '\data_exp_2\cefixime\*.mat']);
load([files.folder '\cefixime_data.mat'])
fittedVals= truncated_fit_simultaneouslyBetaLactams(NGVals,antibioticConcVals,'cefixime',phi_max(6),-1.5);
MICest(6)=fittedVals(1);
Hill_k(6)=fittedVals(2);
pi_min(6)=fittedVals(3);

files = dir([current_dir '\data_exp_2\chloramphenicol\*.mat']);
load([files.folder '\chloramphenicol_data.mat'])
fittedVals= truncated_fit_simultaneouslyChloramphenicol(NGVals,antibioticConcVals,'chloramphenicol',phi_max(7),-0.5);
MICest(7)=fittedVals(1);
Hill_k(7)=fittedVals(2);
pi_min(7)=fittedVals(3);

files = dir([current_dir '\data_exp_2\tetracycline\*.mat']);
load([files.folder '\tetracycline_data.mat'])
fittedVals= truncated_fit_simultaneously(NGVals,antibioticConcVals,'tetracycline',phi_max(8),-0.0005);
MICest(8)=fittedVals(1);
Hill_k(8)=fittedVals(2);
pi_min(8)=fittedVals(3);


files = dir([current_dir '\data_exp_2\ciprofloxacin\*.mat']);
load([files.folder '\ciprofloxacin_data.mat'])
fittedVals= truncated_fit_simultaneouslyCipro(NGVals,antibioticConcVals,'ciprofloxacin',phi_max(9),-1);
MICest(9)=fittedVals(1);
Hill_k(9)=fittedVals(2);
pi_min(9)=fittedVals(3);

results_exp_2 = [phi_max,pi_min,MICest,Hill_k];