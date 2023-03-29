clear
close all
load('500mg_times_6_8h_apart')
colourVals={'b','r',[255,140,0]./255};
lineStyle={'-','--','-.'};
testingMIC_index=5;
cut_off_AUC=find(areaAboveMIC_in(testingMIC_index,:)>147 & areaAboveMIC_in(testingMIC_index,:)<150);%samples that had cut-off MIC near 148-150. 
clearedSamples=cut_off_AUC(find(clearence_Times(cut_off_AUC,testingMIC_index)<7));%AUC/MIC index near the cut-off of 150 but cleared
unclearedSamples=cut_off_AUC(find(clearence_Times(cut_off_AUC,testingMIC_index)>7));%AUC/MIC index near the cut-off of 150 but did not cleared


for i=1:3

plot(timeRunVal{testingMIC_index,unclearedSamples(i)}/24,total_BacteriaVal{testingMIC_index,unclearedSamples(i)},'color',colourVals{i},'linestyle',lineStyle{i},'linewidth',2);%total bacterial load over time using MyPara
ylabel('total NG (bacteria)')
xlabel('time (days)')
hold on

set(gca,'yscale','log','fontsize',14,'Box','off')
areaAboveMIC_in(testingMIC_index,unclearedSamples(i))


end

legend({'147','148','149'},'Box','off');
plot(0:45,10*ones(length(0:45),1),'k','linestyle',':','linewidth',2,'Handlevisibility','off')
ylim([5 5e7])
xlim([0 45])

x0=10;
y0=10;
width=500;
height=300;
set(gcf,'units','points','position',[x0,y0,width,height])

annotation('textbox',...
    [0.0979162995594713 0.9025 0.0283685756240823 0.0874999999999998],...
    'String',{'(a)'},...
    'LineStyle','none',...
    'FontSize',14,...
    'FitBoxToText','off');




