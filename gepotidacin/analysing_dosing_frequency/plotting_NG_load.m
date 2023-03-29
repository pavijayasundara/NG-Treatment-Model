function plotting_NG_load(total_BacteriaVal,colour)

qunatileVals_total_NG=NaN(size(total_BacteriaVal,1),5);

for i=1:size(total_BacteriaVal,1)
    a=total_BacteriaVal(i,:);
    qunatileVals_total_NG(i,:)=quantile(a,[0.25,0.5,0.75,0.025,0.975]);
end
qunatileVals_total_NG(sum(isnan(qunatileVals_total_NG), 2) == 5, :) = [];%

start=1;
stop=145;
linewidth=3;
numRowsPlotting=start:stop;

plot(numRowsPlotting/24,qunatileVals_total_NG(start:stop,2),'LineWidth',linewidth,'color',colour);
 set(gca,'fontsize',18,'Box','off','yscale','log')
xlabel('time (days)')
ylabel('total NG (bacteria)')
set(gca,'fontsize',14)

end



