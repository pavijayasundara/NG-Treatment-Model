function phi_max=phi_max_all_antibioitcs()
clear
close all
numAntibioitcs=8;
phi_max=NaN(numAntibioitcs,1);
current_dir = pwd;


files = dir([current_dir '\data_exp_1\gentamicin\*.mat']);
load([files.folder '\gentamicin_data.mat'])
gentamicin_phiMax= fitToFoerster_phiMax(NGVals,'gentamicin');
phi_max(1)=gentamicin_phiMax;

files = dir([current_dir '\data_exp_1\spectinomycin\*.mat']);
load([files.folder '\spec_data.mat'])

spectinomycin_phiMax= fitToFoerster_phiMax(NGVals,'spectinomycin');
phi_max(2)=spectinomycin_phiMax;

files = dir([current_dir '\data_exp_1\azithromycin\*.mat']);
load([files.folder '\azithro_data.mat'])

azithromycin_phiMax= fitToFoerster_phiMax(NGVals,'azithromycin');
phi_max(3)=azithromycin_phiMax;

files = dir([current_dir '\data_exp_1\penicillin\*.mat']);
load([files.folder '\penicillin_data.mat'])
penicillin_phiMax= fitToFoerster_phiMax(NGVals,'penicillin');
phi_max(4)=penicillin_phiMax;

files = dir([current_dir '\data_exp_1\ceftriaxone\*.mat']);
load([files.folder '\ceftriaxone_data.mat'])
ceftriaxone_phiMax= fitToFoerster_phiMax(NGVals,'ceftriaxone');
phi_max(5)=ceftriaxone_phiMax;

files = dir([current_dir '\data_exp_1\cefixime\*.mat']);
load([files.folder '\cefixime_data.mat'])
cefixime_phiMax= fitToFoerster_phiMax(NGVals,'cefixime');
phi_max(6)=cefixime_phiMax;

files = dir([current_dir '\data_exp_1\chloramphenicol\*.mat']);
load([files.folder '\chloramphenicol_data.mat'])
chloramphenicol_phiMax= fitToFoerster_phiMax(NGVals,'chloramphenicol');
phi_max(7)=chloramphenicol_phiMax;

files = dir([current_dir '\data_exp_1\tetracycline\*.mat']);
load([files.folder '\tetracycline_data.mat'])

tetracycline_phiMax= fitToFoerster_phiMax(NGVals,'tetracycline');
phi_max(8)=tetracycline_phiMax;


files = dir([current_dir '\data_exp_1\ciprofloxacin\*.mat']);
load([files.folder '\ciprofloxacin_data.mat'])
ciprofloxacin_phiMax= fitToFoerster_phiMax(NGVals,'ciprofloxacin');
phi_max(9)=ciprofloxacin_phiMax;

end