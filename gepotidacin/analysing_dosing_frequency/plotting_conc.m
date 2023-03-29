function [h2,qunatileVals_antibioitc]=plotting_conc(antibiotic,colour,adjust)
linewidth=3;
qunatileVals_antibioitc=NaN(size(antibiotic,1),5);

for i=1:size(antibiotic,1)
    qunatileVals_antibioitc(i,:)=quantile(antibiotic(i,:),[0.25,0.5,0.75,0.025,0.975]);
end
qunatileVals_antibioitc(sum(isnan(qunatileVals_antibioitc), 2) == 5, :) = [];
qunatileVals_antibioitc= qunatileVals_antibioitc(any(qunatileVals_antibioitc,2),:);
start=1;
stop=find(qunatileVals_antibioitc(:,4)<0,1)-1;
numRowsPlotting=start+88:stop+88;
 numRowsPlotting=numRowsPlotting-adjust;
h2=plot((numRowsPlotting/24),qunatileVals_antibioitc(start:stop,2),'LineWidth',linewidth,'color',colour);
 set(gca,'fontsize',18,'Box','off','yscale','log')

end