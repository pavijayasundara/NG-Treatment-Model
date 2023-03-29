
clear
close all
fs=20;%font size

figure('units','normalized','outerposition',[0 0 1 1])
testingMIC_index=1;
load('240_genta_1g_AZ_conc')

clearvars areaAboveMIC_ex areaAboveMIC_in 

[areaAboveMIC_ex]=tot_area_under_curve_multiple_doses(Loewe_ExtraVal,potent_ex_val);
[areaAboveMIC_in]=tot_area_under_curve_multiple_doses(Loewe_IntraVal,potent_in_val);

[binAUC_in_clear,binAUC_in_unclear]=numberAUC_bin(areaAboveMIC_in,clearence_Times,40,60,20,120,testingMIC_index);%%bin start, end, increment,final Val
[binAUC_ex_clear,binAUC_ex_unclear]=numberAUC_bin(areaAboveMIC_ex,clearence_Times,40,60,20,120,testingMIC_index);%%bin start, end, increment,final Val

AUCBar_in(:,1)=binAUC_in_clear;
AUCBar_in(:,2)=binAUC_in_unclear;
AUCBar_ex(:,1)=binAUC_ex_clear;
AUCBar_ex(:,2)=binAUC_ex_unclear;


load('80times_3_8hGEN_1g_AZ')
clearvars areaAboveMIC_ex_gen areaAboveMIC_in_gen timeAboveMIC_in timeAboveMIC_ex

[areaAboveMIC_ex]=tot_area_under_curve_multiple_doses(Loewe_ExtraVal,potent_ex_val);
[areaAboveMIC_in]=tot_area_under_curve_multiple_doses(Loewe_IntraVal,potent_in_val);

[binAUC_in_clear,binAUC_in_unclear]=numberAUC_bin(areaAboveMIC_in,clearence_Times,40,60,20,120,testingMIC_index);%%bin start, end, increment,final Val
[binAUC_ex_clear,binAUC_ex_unclear]=numberAUC_bin(areaAboveMIC_ex,clearence_Times,40,60,20,120,testingMIC_index);%%bin start, end, increment,final Val

AUCBar_in(:,3)=binAUC_in_clear;
AUCBar_in(:,4)=binAUC_in_unclear;
AUCBar_ex(:,3)=binAUC_ex_clear;
AUCBar_ex(:,4)=binAUC_ex_unclear;

load('120times_2_8hGEN_1g_AZ')
clearvars areaAboveMIC_ex_gen areaAboveMIC_in_gen timeAboveMIC_in timeAboveMIC_ex
[areaAboveMIC_ex]=tot_area_under_curve_multiple_doses(Loewe_ExtraVal,potent_ex_val);
[areaAboveMIC_in]=tot_area_under_curve_multiple_doses(Loewe_IntraVal,potent_in_val);

[binAUC_in_clear,binAUC_in_unclear]=numberAUC_bin(areaAboveMIC_in,clearence_Times,40,60,20,120,testingMIC_index);%%bin start, end, increment,final Val
[binAUC_ex_clear,binAUC_ex_unclear]=numberAUC_bin(areaAboveMIC_ex,clearence_Times,40,60,20,120,testingMIC_index);%%bin start, end, increment,final Val

AUCBar_in(:,5)=binAUC_in_clear;
AUCBar_in(:,6)=binAUC_in_unclear;
AUCBar_ex(:,5)=binAUC_ex_clear;
AUCBar_ex(:,6)=binAUC_ex_unclear;

groupLabels = categorical({'60-80' '80-100' '100-120','120-140','140-160'}); 

s1=subplot(1,2,1);
h=bar(groupLabels,AUCBar_in);
colouBar(h)
xlabel('AUC/MIC_{in} (h)')
ylim([0 5500])
set(gca,'fontsize',fs,'yscale','log','ytick',[10,100,1000],...
    'yticklabel',{'10^1','10^2','10^3'},'Box', 'off')
groupLabels = categorical({'60-80' '80-100' '100-120','120-140','140-160'}); 

subplot(1,2,2)  
b=bar(groupLabels,AUCBar_ex);
colouBar(b)
xlabel('AUC/MIC_{ex} (h)')
ylim([0 5500])
set(gca,'fontsize',fs,'yscale','log','ytick',[1,10,100,1000],...
    'yticklabel',{'10^0','10^1','10^2','10^3'},'Box', 'off')

legend([h(1,1),h(1,2),h(1,3),h(1,4),h(1,5),h(1,6)],{'500mg×6, 8h uncleared','500mg×6, 8h cleared'...
    '1500mgmg×2, 24h uncleared','1500mgmg×2, 24h cleared','500mg×6, 12h uncleared'...
    '500mg×6, 12h cleared'},... 
    'NumColumns',3,...
    'EdgeColor',[1 1 1],...
    'Position',[0.254636415847194 0.0181771679163832 0.496339666306711 0.0237173556219415],...
    'Box','off')

print(gcf, 'test_new', '-dpng', '-r310' ); 