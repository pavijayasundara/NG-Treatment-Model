function [fittedVals,errorVals,confInt,h3]= fitting_mean_effect2(Hilleffect,antibioticConcVals,titleName,phi_minIntGuess)

conc=antibioticConcVals;

Y=Hilleffect; % data on the mean hill function effect

modelHill=@(x,conc)x(1)-(((x(1)-x(4)).*(conc./x(2)).^x(3))./((conc./x(2)).^x(3)-(x(4)/x(1))));

%parameters estimated: phi-max,MIC and Hill_k,phi_min
initial=[Y(1);0.004;0.2;phi_minIntGuess];%initial guesses of the estimates
lb=[0 0 0 -10];%lower bound of the estimates
ub=[1 10 10 0];%upper bound of the estimates

%options=optimoptions(@lsqnonlin,'MaxFunctionEvaluations',10000,'MaxIterations',10000);
[fittedVals,errorVals,residual,~,~,~,jacobian]= lsqcurvefit(modelHill,initial,conc,Y,lb,ub);

fittedValues=modelHill(fittedVals',conc);%columns are time points, rows are concentrations 

confInt=nlparci(fittedVals,residual,'jacobian',jacobian);

scatter(conc,Y,20,'k','filled');
hold on
h3=plot(conc,fittedValues,'k--');
        

hold off

%legend(plotHandles, plotLabels,'Location','NorthEastOutside','Box','off');
end

