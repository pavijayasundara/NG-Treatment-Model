function parameterEstimates= fitToFoerster_phiMax(NGVals,titleName)
%fitting to NG data of Foerster. Estimating phi_max by using the
%concentration without antibiotic.
%close all

    parameterEstimates=NaN(1,1);%12 concentrations and three parameter estimates.
    ErrorEstimates=NaN(1,1);


    initial=0.5;

    thours=[0;1;2;3;4;5;6];

    %-4 hour time point is not considered
    %data=log10(NGVals(j,2:end));%NG data in log
    data=NGVals(1,2:end);%NG data in log
    data=data';
    lb=0;
    ub=1;

    x=[thours,data];
    [xmulti,errormulti]=fmincon(@(parameters)ssq_phiMax(parameters,x),initial,[],[],[],[],lb,ub,[]);
    fittedValues=modelphiMax(xmulti',thours,NGVals(1,2));
    parameterEstimates(1,:)= xmulti;
    ErrorEstimates(1)=errormulti;
   
    scatter(thours,data,20,'or','filled')
    hold on

    plot(thours,fittedValues,'b')
    title(titleName)
    ylabel('NG')
    xlabel('time(h)')
    set(gca,'yScale','log','fontsize',14)
    
end

 


