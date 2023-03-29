clear
close all
%load('FoersterEstimatesofParavalues.mat');
%load('MYestimatesusingMLE.mat');
totSSE=NaN(9,1);
SSE_conc=NaN(11,18);
x1=3;
x2=3;
figure('units','normalized','outerposition',[0 0 1 1])
names=["gentamicin","spectinomycin","azithromycin","penicillin","ceftriaxone","cefixime","chloramphenicol","tetracycline","ciprofloxain"];

%my estimates
 Hill_kest=[1.702922562;0.997854233;0.910375802;1.095883648;1.751173383;1.205518482;1.080105605;0.87024053;1.001038367];
 phi_maxest=[0.86716845;0.7516419;0.683800635;0.723164057;0.713691958;0.735962308;0.676470001;0.761810405;0.753514106];
 phi_minest=[-8.179497977;-9.457113037;-1.498685512;-1.277131048;-0.451975325;-0.507544298;-0.245195858;-0.184498809;-9.999999999];
 MICest=[0.236274166;7.423857592;0.026728278;0.003430895;0.000296956;0.000212656;0.55394414;0.567921636;0.001929046];


%Foesrter study estimates
% Hill_kest=[0.820000000000000;2.41000000000000;2.43000000000000;1.19000000000000;1.69000000000000;2.07000000000000;1.54000000000000;1.14000000000000;1.04];
% phi_maxest=[0.860000000000000;0.760000000000000;0.670000000000000;0.770000000000000;0.700000000000000;0.890000000000000;0.850000000000000;0.720000000000000;0.98];
% phi_minest=[-206.800000000000;-10.3000000000000;-2.25000000000000;-2.06000000000000;-0.740000000000000;-0.640000000000000;-0.120000000000000;-0.250000000000000;-7.34];
% MICest=[0.1522;5.64520000000000;0.0267000000000000;0.00530000000000000;0.000200000000000000;0.000100000000000000;0.376700000000000;0.325900000000000;0.0017];

%Foerster other experiment data
% Hill_kest=[1.17;1.61;2.64;1.01;1.58;1.42;2.04;0.9;1.19];
% phi_maxest=[0.96;0.71;0.6;1.05;0.8;0.75;0.61;0.83;0.43];
% phi_minest=[-7.96;-8.94;-2.07;-1.15;-0.46;-0.87;-0.1;-0.13;-10.4];
% MICest=[0.2117;4.6836;0.0238;0.0029;0.0004;0.0004;0.5762;0.7067;0.0018];

%Foerster mean parameter values in the main text
% Hill_kest=[1;2;2.5;1.1;1.6;1.7;1.8;1;1.1];
% phi_maxest=[0.9	;0.7;0.6;0.9;0.8	;0.8;0.7;0.8;0.7];
% phi_minest=[-106.9;-9.6;	-2.2;	-1.6	;-0.6	;-0.8	;-0.1	;-0.2;	-8.9];
% MICest=[0.2	;5	;0.03	;0.004;	0.0003;	0.0002;	0.5	;0.5;	0.002];

for i=1:9
subplot(x1,x2,i)
if i==1
load('C:\Users\z5137734\OneDrive - UNSW\PhD\treatment model\treatment model with new equations\Foerester model reproduce\fittingMyModelToFoerster\data_3\gentamicin\gentamicin_data')
end
if i==2
load('C:\Users\z5137734\OneDrive - UNSW\PhD\treatment model\treatment model with new equations\Foerester model reproduce\fittingMyModelToFoerster\data_3\spec\spec_data')
end
if i==3
load('C:\Users\z5137734\OneDrive - UNSW\PhD\treatment model\treatment model with new equations\Foerester model reproduce\fittingMyModelToFoerster\data_3\azithromycin\azithroData')
end
if i==4
load('C:\Users\z5137734\OneDrive - UNSW\PhD\treatment model\treatment model with new equations\Foerester model reproduce\fittingMyModelToFoerster\data_3\penicillin\penicillin_data')
end
if i==5
    load('C:\Users\z5137734\OneDrive - UNSW\PhD\treatment model\treatment model with new equations\Foerester model reproduce\fittingMyModelToFoerster\data_3\ceftriaxone\ceftriaxone_data')
end
if i==6
load('C:\Users\z5137734\OneDrive - UNSW\PhD\treatment model\treatment model with new equations\Foerester model reproduce\fittingMyModelToFoerster\data_3\cefixime\cefixime_data')
end
if i==7
    load('C:\Users\z5137734\OneDrive - UNSW\PhD\treatment model\treatment model with new equations\Foerester model reproduce\fittingMyModelToFoerster\data_3\chloro\chloramphenicol_data')
end
if i==8
    load('C:\Users\z5137734\OneDrive - UNSW\PhD\treatment model\treatment model with new equations\Foerester model reproduce\fittingMyModelToFoerster\data_3\tetra\tetracycline_data')
end

if i==9
    load('C:\Users\z5137734\OneDrive - UNSW\PhD\treatment model\treatment model with new equations\Foerester model reproduce\fittingMyModelToFoerster\data_3\cipro\ciprofloxacin_data')
end

if i==1
[SSEconc,totalSSE]= SSEFoesrsterGenta(NGVals,antibioticConcVals,names(i),phi_maxest(i),phi_minest(i),MICest(i),Hill_kest(i));
SSEconc(9:11)=NaN;
elseif i==2
[SSEconc,totalSSE]= SSEFoesrsterGenta(NGVals,antibioticConcVals,names(i),phi_maxest(i),phi_minest(i),MICest(i),Hill_kest(i));
SSEconc(9:11)=NaN;    
elseif i==4
[SSEconc,totalSSE]= SSEFoesrsterPenicillin(NGVals,antibioticConcVals,names(i),phi_maxest(i),phi_minest(i),MICest(i),Hill_kest(i));
SSEconc(9:11)=NaN;      
elseif i==5
[SSEconc,totalSSE]= SSEFoesrsterBetaLactams(NGVals,antibioticConcVals,names(i),phi_maxest(i),phi_minest(i),MICest(i),Hill_kest(i));
SSEconc(9:11)=NaN;
elseif i==6
[SSEconc,totalSSE]= SSEFoesrsterBetaLactams(NGVals,antibioticConcVals,names(i),phi_maxest(i),phi_minest(i),MICest(i),Hill_kest(i));
SSEconc(9:11)=NaN;
elseif i==7
[SSEconc,totalSSE]= SSEFoesrsterChloramphenicol(NGVals,antibioticConcVals,names(i),phi_maxest(i),phi_minest(i),MICest(i),Hill_kest(i));
SSEconc(11)=NaN;    
elseif i==9
[SSEconc,totalSSE]= SSEFoesrsterCipro(NGVals,antibioticConcVals,names(i),phi_maxest(i),phi_minest(i),MICest(i),Hill_kest(i));
SSEconc(11)=NaN;       
else

[SSEconc,totalSSE]= SSEFoesrsterEstimates2(NGVals,antibioticConcVals,names(i),phi_maxest(i),phi_minest(i),MICest(i),Hill_kest(i));

end
% SSE_conc(:,2*i-1)=antibioticConcVals(2:end);
% SSE_conc(:,2*i)=SSEconc;

totSSE(i,1)=totalSSE;

end

% subplot(x1,x2,2)
% %clear
% 
% [totalSSE,SSEconc]= SSEFoesrsterEstimates(NGVals,antibioticConcVals,'spectinomycin',phi_max(12),phi_min(2),MIC(2),Hill_k(2));
% SSE_conc2=NaN(11,2);%errors at each concnetration. 11 concentrations 
% SSE_conc2(:,1)=antibioticConcVals(2:end);
% SSE_conc2(:,2)=SSEconc;
% totSSE(2)=totalSSE;
% 
% 
%