function phi_max=phi_max_all_antibioitcs()
clear
close all
numAntibioitcs=8;
phi_max=NaN(numAntibioitcs,1);
figure('units','normalized','outerposition',[0 0 1 1])

subplot(2,4,1)
load('C:\Users\z5137734\OneDrive - UNSW\PhD\treatment model\treatment model with new equations\Foerester model reproduce\fittingMyModelToFoerster\data_3\gentamicin\gentamicin_data')
gentamicin_phiMax= fitToFoerster_phiMax(NGVals,'gentamicin');
phi_max(1)=gentamicin_phiMax;

subplot(2,4,2)
%clear
load('C:\Users\z5137734\OneDrive - UNSW\PhD\treatment model\treatment model with new equations\Foerester model reproduce\fittingMyModelToFoerster\data_3\spec\spec_data')
spectinomycin_phiMax= fitToFoerster_phiMax(NGVals,'spectinomycin');
phi_max(2)=spectinomycin_phiMax;

subplot(2,4,3)
%clear
load('C:\Users\z5137734\OneDrive - UNSW\PhD\treatment model\treatment model with new equations\Foerester model reproduce\fittingMyModelToFoerster\data_3\azithromycin\azithroData')
azithromycin_phiMax= fitToFoerster_phiMax(NGVals,'azithromycin');
phi_max(3)=azithromycin_phiMax;

subplot(2,4,4)
%clear
load('C:\Users\z5137734\OneDrive - UNSW\PhD\treatment model\treatment model with new equations\Foerester model reproduce\fittingMyModelToFoerster\data_3\penicillin\penicillin_data')
penicillin_phiMax= fitToFoerster_phiMax(NGVals,'penicillin');
phi_max(4)=penicillin_phiMax;

subplot(2,4,5)
%clear
load('C:\Users\z5137734\OneDrive - UNSW\PhD\treatment model\treatment model with new equations\Foerester model reproduce\fittingMyModelToFoerster\data_3\ceftriaxone\ceftriaxone_data')
ceftriaxone_phiMax= fitToFoerster_phiMax(NGVals,'ceftriaxone');
phi_max(5)=ceftriaxone_phiMax;

subplot(2,4,6)
%clear
load('C:\Users\z5137734\OneDrive - UNSW\PhD\treatment model\treatment model with new equations\Foerester model reproduce\fittingMyModelToFoerster\data_3\cefixime\cefixime_data')
cefixime_phiMax= fitToFoerster_phiMax(NGVals,'cefixime');
phi_max(6)=cefixime_phiMax;

subplot(2,4,7)
%clear
load('C:\Users\z5137734\OneDrive - UNSW\PhD\treatment model\treatment model with new equations\Foerester model reproduce\fittingMyModelToFoerster\data_3\chloro\chloramphenicol_data')
chloramphenicol_phiMax= fitToFoerster_phiMax(NGVals,'chloramphenicol');
phi_max(7)=chloramphenicol_phiMax;

subplot(2,4,8)
%clear
load('C:\Users\z5137734\OneDrive - UNSW\PhD\treatment model\treatment model with new equations\Foerester model reproduce\fittingMyModelToFoerster\data_3\tetra\tetracycline_data')
tetracycline_phiMax= fitToFoerster_phiMax(NGVals,'tetracycline');
phi_max(8)=tetracycline_phiMax;

figure

%clear
load('C:\Users\z5137734\OneDrive - UNSW\PhD\treatment model\treatment model with new equations\Foerester model reproduce\fittingMyModelToFoerster\data_3\cipro\ciprofloxacin_data')
ciprofloxacin_phiMax= fitToFoerster_phiMax(NGVals,'ciprofloxacin');
phi_max(9)=ciprofloxacin_phiMax;

end