%%analysing dosing frequency%%%

%%%finding the indexes that cleared in < 7 days in 8 and 12h apart but
%%%failed% at 24 h

testingMIC_index=5;
close all

load('frequency data')
indexVals=find(c1<7 & c2<7 & c3>7);
indexVals=indexVals(2);
x0=10;
y0=10;
width=1050;
height=310;
set(gcf,'units','points','position',[x0,y0,width,height])

subplot(1,2,1)
[h1,~]=plotting_conc(a1_in(88:end,indexVals),'r',0);
hold on
[h2,~]=plotting_conc(a2_in(88:end,indexVals),'b',0);
hold on
[h3,~]=plotting_conc(a3_in(88:end,indexVals),[255,140,0]./255,0);

xlim([0 10])
ylabel(['intracellular       ';'concentration (mg/L)']);
xlabel('time(days)')
set(gca,'fontsize',18,'ytick',[0.5 1 2 3 4 5],'xtick',[0,1,2,3,4,5,6,7,8,9,10])
ylim([0.25 6])
legend([h1,h2,h3],{'500 × 6, 8h','500 × 6, 12h','500 × 6, 24h'},'Box','off');

subplot(1,2,2)
plotting_NG_load(b1(:,indexVals),'r')
hold on
plotting_NG_load(b2(:,indexVals),'b')
hold on 
plot(timeRunVal(:,indexVals)./24,b3(:,indexVals),'linewidth',3,'color',[255,140,0]./255)
hold on
plot(0:12,10*ones(length(0:12),1),'k','linestyle','--','linewidth',3)
ylabel('total NG (bacteria)')
xlabel('time (days)')
set(gca,'fontsize',18,'ytick',[10 10^2 10^3 10^4 10^5 10^6 10^7],'xtick',[0,1,2,3,4,5,6,7,8,9,10])
xlim([0 10])
ylim([1 10^8])

annotation('textbox',...
    [0.101 0.934814814814815 0.0182592592592593 0.037037037037037],...
    'String',{'(a)'},...
    'LineStyle','none',...
    'FontSize',18,...
    'FitBoxToText','off');

annotation('textbox',...
    [0.540259259259259 0.927407407407407 0.0182592592592592 0.0370370370370371],...
    'String','(b)',...
    'LineStyle','none',...
    'FontSize',18,...
    'FitBoxToText','off');
annotation('line',[0.130690161527166 0.462555066079295],...
    [0.503830917874396 0.502415458937198],'LineWidth',3,'LineStyle','--');


