clear
close all
fs=100;%font size

figure
set(gcf, 'PaperPosition', [0 0 195 180])

%%testing the cut-off required for MIC=1 mg/L
testingMIC_index=5;
load('500mg_times_6_8h_apart')
clearvars areaAboveMIC_ex areaAboveMIC_in timeAboveMIC_in timeAboveMIC_ex

[areaAboveMIC_ex,~]=tot_area_under_curve_multiple_doses(antibiotic,timeRun,6,8,LHS_treatment);
[areaAboveMIC_in,~]=tot_area_under_curve_multiple_doses(intracellular_atibioitc,timeRun,6,8,LHS_treatment);
[timeAboveMIC_in,~]=time_above_MIC_multiple_doses(intracellular_atibioitc,timeRun,6,8);
[timeAboveMIC_ex,~]=time_above_MIC_multiple_doses(antibiotic,timeRun,6,8);

[binAUC_in_clear,binAUC_in_unclear]=numberAUC_bin(areaAboveMIC_in,clearence_Times,110,130,20,170);%%bin start, end, increment,final Val
[binAUC_ex_clear,binAUC_ex_unclear]=numberAUC_bin(areaAboveMIC_ex,clearence_Times,70,90,20,130);%%bin start, end, increment,final Val
[bin_tMIC_in_clear,bin_tMIC_in_unclear]=numberAUC_bin(timeAboveMIC_in,clearence_Times,35,45,10,65);%%bin start, end, increment,final Val
[bin_tMIC_ex_clear,bin_tMIC_ex_unclear]=numberAUC_bin(timeAboveMIC_ex,clearence_Times,35,45,10,65);%%bin start, end, increment,final Val

[CmaxAboveMIC_in,~]=Cmax(intracellular_atibioitc,timeRun);
[CmaxAboveMIC_ex,~]=Cmax(antibiotic,timeRun);
[bin_Cmax_in_clear,bin_Cmax_in_unclear]=numberAUC_bin(CmaxAboveMIC_in,clearence_Times,1,3,2,7);%%bin start, end, increment,final Val
[bin_Cmax_ex_clear,bin_Cmax_ex_unclear]=numberAUC_bin(CmaxAboveMIC_ex,clearence_Times,1,3,2,7);%%bin start, end, increment,final Val

AUCBar_in(:,1)=binAUC_in_unclear;
AUCBar_in(:,2)=binAUC_in_clear;
AUCBar_ex(:,1)=binAUC_ex_unclear;
AUCBar_ex(:,2)=binAUC_ex_clear;
t_MIC_Bar_in(:,1)=bin_tMIC_in_unclear;
t_MIC_Bar_in(:,2)=bin_tMIC_in_clear;
t_MIC_Bar_ex(:,1)=bin_tMIC_ex_unclear;
t_MIC_Bar_ex(:,2)=bin_tMIC_ex_clear;
C_max_Bar_in(:,1)=bin_Cmax_in_unclear;
C_max_Bar_in(:,2)=bin_Cmax_in_clear;
C_max_Bar_ex(:,1)=bin_Cmax_ex_unclear;
C_max_Bar_ex(:,2)=bin_Cmax_ex_clear;

load('1500mg_times_2_24h_apart')
clearvars areaAboveMIC_ex areaAboveMIC_in timeAboveMIC_in timeAboveMIC_ex

[areaAboveMIC_ex,~]=tot_area_under_curve_multiple_doses(antibiotic,timeRun,2,24,LHS_treatment);
[areaAboveMIC_in,~]=tot_area_under_curve_multiple_doses(intracellular_atibioitc,timeRun,2,24,LHS_treatment);
[timeAboveMIC_in,~]=time_above_MIC_multiple_doses(intracellular_atibioitc,timeRun,2,24);
[timeAboveMIC_ex,~]=time_above_MIC_multiple_doses(antibiotic,timeRun,2,24);

[binAUC_in_clear,binAUC_in_unclear]=numberAUC_bin(areaAboveMIC_in,clearence_Times,110,130,20,170);%%bin start, end, increment,final Val
[binAUC_ex_clear,binAUC_ex_unclear]=numberAUC_bin(areaAboveMIC_ex,clearence_Times,70,90,20,130);%%bin start, end, increment,final Val
[bin_tMIC_in_clear,bin_tMIC_in_unclear]=numberAUC_bin(timeAboveMIC_in,clearence_Times,35,45,10,65);%%bin start, end, increment,final Val
[bin_tMIC_ex_clear,bin_tMIC_ex_unclear]=numberAUC_bin(timeAboveMIC_ex,clearence_Times,35,45,10,65);%%bin start, end, increment,final Val

[CmaxAboveMIC_in,~]=Cmax(intracellular_atibioitc,timeRun);
[CmaxAboveMIC_ex,~]=Cmax(antibiotic,timeRun);
[bin_Cmax_in_clear,bin_Cmax_in_unclear]=numberAUC_bin(CmaxAboveMIC_in,clearence_Times,1,3,2,7);%%bin start, end, increment,final Val
[bin_Cmax_ex_clear,bin_Cmax_ex_unclear]=numberAUC_bin(CmaxAboveMIC_ex,clearence_Times,1,3,2,7);%%bin start, end, increment,final Val

AUCBar_in(:,3)=binAUC_in_unclear;
AUCBar_in(:,4)=binAUC_in_clear;
AUCBar_ex(:,3)=binAUC_ex_unclear;
AUCBar_ex(:,4)=binAUC_ex_clear;
t_MIC_Bar_in(:,3)=bin_tMIC_in_unclear;
t_MIC_Bar_in(:,4)=bin_tMIC_in_clear;
t_MIC_Bar_ex(:,3)=bin_tMIC_ex_unclear;
t_MIC_Bar_ex(:,4)=bin_tMIC_ex_clear;
C_max_Bar_in(:,3)=bin_Cmax_in_unclear;
C_max_Bar_in(:,4)=bin_Cmax_in_clear;
C_max_Bar_ex(:,3)=bin_Cmax_ex_unclear;
C_max_Bar_ex(:,4)=bin_Cmax_ex_clear;


load('500mg_times_6_12h_apart')
clearvars areaAboveMIC_ex areaAboveMIC_in timeAboveMIC_in timeAboveMIC_ex
[areaAboveMIC_ex,summary_area_above_MIC_ex]=tot_area_under_curve_multiple_doses(antibiotic,timeRun,6,12,LHS_treatment);
[areaAboveMIC_in,summary_area_above_MIC_in]=tot_area_under_curve_multiple_doses(intracellular_atibioitc,timeRun,6,12,LHS_treatment);
[timeAboveMIC_in,summary_t_above_MIC_in]=time_above_MIC_multiple_doses(intracellular_atibioitc,timeRun,6,12);
[timeAboveMIC_ex,summary_t_above_MIC_ex]=time_above_MIC_multiple_doses(antibiotic,timeRun,6,12);


[binAUC_in_clear,binAUC_in_unclear]=numberAUC_bin(areaAboveMIC_in,clearence_Times,110,130,20,170);%%bin start, end, increment,final Val
[binAUC_ex_clear,binAUC_ex_unclear]=numberAUC_bin(areaAboveMIC_ex,clearence_Times,70,90,20,130);%%bin start, end, increment,final Val
[bin_tMIC_in_clear,bin_tMIC_in_unclear]=numberAUC_bin(timeAboveMIC_in,clearence_Times,35,45,10,65);%%bin start, end, increment,final Val
[bin_tMIC_ex_clear,bin_tMIC_ex_unclear]=numberAUC_bin(timeAboveMIC_ex,clearence_Times,35,45,10,65);%%bin start, end, increment,final Val

[CmaxAboveMIC_in,~]=Cmax(intracellular_atibioitc,timeRun);
[CmaxAboveMIC_ex,~]=Cmax(antibiotic,timeRun);
[bin_Cmax_in_clear,bin_Cmax_in_unclear]=numberAUC_bin(CmaxAboveMIC_in,clearence_Times,1,3,2,7);%%bin start, end, increment,final Val
[bin_Cmax_ex_clear,bin_Cmax_ex_unclear]=numberAUC_bin(CmaxAboveMIC_ex,clearence_Times,1,3,2,7);%%bin start, end, increment,final Val

AUCBar_in(:,5)=binAUC_in_unclear;
AUCBar_in(:,6)=binAUC_in_clear;
AUCBar_ex(:,5)=binAUC_ex_unclear;
AUCBar_ex(:,6)=binAUC_ex_clear;
t_MIC_Bar_in(:,5)=bin_tMIC_in_unclear;
t_MIC_Bar_in(:,6)=bin_tMIC_in_clear;
t_MIC_Bar_ex(:,5)=bin_tMIC_ex_unclear;
t_MIC_Bar_ex(:,6)=bin_tMIC_ex_clear;
C_max_Bar_in(:,5)=bin_Cmax_in_unclear;
C_max_Bar_in(:,6)=bin_Cmax_in_clear;
C_max_Bar_ex(:,5)=bin_Cmax_ex_unclear;
C_max_Bar_ex(:,6)=bin_Cmax_ex_clear;



groupLabels = categorical({'110-130' '130-150' '150-170','170-190'});
s1=subplot(3,2,1);
h=bar(groupLabels,AUCBar_in);
colouBar(h)
xlabel('AUC/MIC_{in} (h)')
ylim([0 5500])
set(gca,'fontsize',fs,'yscale','log','ytick',[1,10,100,1000],...
    'yticklabel',{'10^0','10^1','10^2','10^3'},'Box', 'off')

groupLabels = categorical({'110-130' '130-150' '150-170','170-190'});
subplot(3,2,2)
b=bar(groupLabels,AUCBar_ex);
colouBar(b)
xlabel('AUC/MIC_{ex} (h)')
ylim([0 5500])
set(gca,'fontsize',fs,'yscale','log','ytick',[1,10,100,1000],...
    'yticklabel',{'10^0','10^1','10^2','10^3'},'Box', 'off')

groupLabels = categorical({'35-45' '45-55','55-65','65-75'});
subplot(3,2,3)
b=bar(groupLabels,t_MIC_Bar_in);
colouBar(b)
xlabel('t_{MICin} (h)')
ylim([0 5200])
set(gca,'fontsize',fs,'yscale','log','ytick',[1,10,100,1000],...
    'yticklabel',{'10^0','10^1','10^2','10^3'},'Box', 'off')

groupLabels = categorical({'35-45' '45-55','55-65','65-75'});
subplot(3,2,4)
b=bar(groupLabels,t_MIC_Bar_ex);
colouBar(b)
xlabel('t_{MICex} (h)')
ylim([0 5500])
set(gca,'fontsize',fs,'yscale','log','ytick',[1,10,100,1000],...
    'yticklabel',{'10^0','10^1','10^2','10^3'},'Box', 'off')


groupLabels = categorical({'1-3' '3-5' '5-7','7-9'});
subplot(3,2,5)
b=bar(groupLabels,C_max_Bar_in);
colouBar(b)
xlabel('C_{max}/MIC_{in}')
ylim([0 5500])
set(gca,'fontsize',fs,'yscale','log','ytick',[1,10,100,1000],...
    'yticklabel',{'10^0','10^1','10^2','10^3'},'Box', 'off')

groupLabels = categorical({'1-3' '3-5' '5-7','7-9'});
subplot(3,2,6)
b=bar(groupLabels,C_max_Bar_ex);
colouBar(b)
xlabel('C_{max}/MIC_{ex}')
ylim([0 5500])
set(gca,'fontsize',fs,'yscale','log','ytick',[1,10,100,1000],...
    'yticklabel',{'10^0','10^1','10^2','10^3'},'Box', 'off')


legend([h(1,1),h(1,2),h(1,3),h(1,4),h(1,5),h(1,6)],{'500×6, 8h uncleared','500×6, 8h cleared'...
    '1500mg×2, 24h uncleared','1500mg×2, 24h cleared','500×6, 12h uncleared'...
    '500×6, 12h cleared'},...
    'NumColumns',3,...
    'EdgeColor',[1 1 1],...
    'Position',[0.254636415847194 0.0181771679163832 0.496339666306711 0.0237173556219415],...
    'Box','off')

annotation('textbox',...
    [0.109888888888889 0.639999999999999 0.0286296296296296 0.0414814814814815],...
    'String','(c)',...
    'LineStyle','none',...
    'FontSize',fs,...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);

annotation('textbox',...
    [0.109888888888889 0.928888888888888 0.0286296296296296 0.0414814814814815],...
    'String',{'(a)'},...
    'LineStyle','none',...
    'FontSize',fs,...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);

%%% textbox
annotation('textbox',...
    [0.55062962962963 0.933333333333333 0.0286296296296296 0.0414814814814815],...
    'String','(b)',...
    'LineStyle','none',...
    'FontSize',fs,...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);

%%% textbox
annotation('textbox',...
    [0.547666666666667 0.634074074074073 0.0286296296296297 0.0414814814814815],...
    'String','(d)',...
    'LineStyle','none',...
    'FontSize',fs,...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);

%%% textbox
annotation('textbox',...
    [0.106925925925926 0.331851851851851 0.0286296296296297 0.0414814814814815],...
    'String','(e)',...
    'LineStyle','none',...
    'FontSize',fs,...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);

%%% textbox
annotation('textbox',...
    [0.547666666666667 0.337777777777777 0.0286296296296297 0.0414814814814815],...
    'String','(f)',...
    'LineStyle','none',...
    'FontSize',fs,...
    'FitBoxToText','off',...
    'EdgeColor',[1 1 1]);

annotation('textarrow',[0.0718518518518519 0.391343795228485],...
    [0.705185185185185 0.474701512022598],'TextRotation',90,...
    'String',{'number of samples'},...
    'LineStyle','none',...
    'HeadStyle','none',...
    'FontSize',fs+65);

