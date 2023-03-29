function [fittedVals,errorVals,confInt,SSEconc]= truncated_fit_simultaneouslyChloramphenicol(NGVals,antibioticConcVals,titleName,phi_max,phi_minIntGuess)

thours=[0;1;2;3;4;5;6];
phi_maxVal=phi_max;
conc=antibioticConcVals(2:end);

X=[phi_maxVal;conc;NGVals(2,2)];%concentrations and phi_max value
%cocnvals should be only the 10 concnentratons excluding the 0 and maximum
%concentrations

Y=NGVals(2:end,2:end);

%finding the values which are=100 and replacing them with NaN so that they
%act as missing values. The time variable is unaffected
if find(Y==100,1)>0
        Y(find(Y==100))=NaN;
       
end

%parameters estimated: MIC and Hill_k,phi_min
initial=[0.004;8;phi_minIntGuess];%initial guesses of the estimates
lb=[0 0 -10];%lower bound of the estimates
ub=[10 5 0];%upper bound of the estimates

options=optimoptions(@lsqnonlin,'MaxFunctionEvaluations',10000,'MaxIterations',10000);
[fittedVals,errorVals,residual,exitflag,output,lambda,jacobian]=lsqnonlin(@(xEstimate)errorFun(xEstimate,Y,thours,X),initial,lb,ub,options);
concVals=round(conc,4);
colourVals={[75,0,13]./255,[176,224,230]./255,[0,191,25]./255,[30,144,255]./255,[100,149,237]./255,...
    [70,130,180]./255,[95,158,160]./255,[123,104,238]./255,[72,61,139]./255,[0,0,255]./255,[0,0,128]./255,'g'};
fittedValues=model(fittedVals',thours,X);%columns are time points, rows are concentrations 
errorVals;

%errors at each of the concentrations
err1=log10(fittedValues)-log10(Y);
err1(isnan(err1))=0;%replacing NaN errors by 0 (the points where data are censored)
err2=(err1).^2;
SSEconc=sum(err2,2);%sum of squared errors at each concentration
%95% confidence intervals 
confInt=  nlparci(fittedVals,residual,'jacobian',jacobian);

plotHandles = zeros(1,length(conc));
plotLabels = cell(1,length(conc));
for k=1:length(conc)
scatter(thours,Y(k,:),20,colourVals{k},'filled')
hold on
plotLabels{k}= num2str(concVals(k));
plotHandles(k)=plot(thours,fittedValues(k,:),'color', colourVals{k});
plot(thours,fittedValues(k,:),'color', colourVals{k})
hold on
end
%legend({'data','fit'},'Box','off')
ylabel('NG')
xlabel('time(hours)')
set(gca,'fontsize',14,'yscale','log')
title(titleName)
legend(plotHandles, plotLabels,'Location','NorthEastOutside','Box','off');
end

function err=errorFun(xEstimate,y_meas,thours,X)
    y_est=model(xEstimate,thours,X);
    err=log10(y_meas)-log10(y_est);
    err=err(:);
    err(isnan(err)) = [];
end
function output= model(xEstimate,thours,X)
    
    intVals=repmat(X(end),10,1);%initial value of the ODE system
    
    options = odeset('AbsTol',1e-13,'RelTol',1e-11);
    
    [~,values] = ode45(@(t,y)Equations(t,y,X(1:end-1)),thours,intVals,options);
 
    function s=Equations(~,y,para)

        phi_max=para(1);
        phi_min=xEstimate(3);
        MIC=xEstimate(1);
        Hill_k=xEstimate(2);
        conc=para(2:end,1);%concentrations
    
        s=zeros(10,1);
        s(1)=phi_max*y(1)-y(1)*(((phi_max-phi_min)*(conc(1)/MIC).^Hill_k)/((conc(1)/MIC).^Hill_k-(phi_min/(phi_max)))); 
        s(2)=phi_max*y(2)-y(2)*(((phi_max-phi_min)*(conc(2)/MIC).^Hill_k)/((conc(2)/MIC).^Hill_k-(phi_min/(phi_max)))); 
   
        s(3)=phi_max*y(3)-y(3)*(((phi_max-phi_min)*(conc(3)/MIC).^Hill_k)/((conc(3)/MIC).^Hill_k-(phi_min/(phi_max)))); 
        
        s(4)=phi_max*y(4)-y(4)*(((phi_max-phi_min)*(conc(4)/MIC).^Hill_k)/((conc(4)/MIC).^Hill_k-(phi_min/(phi_max))));
        
        s(5)=phi_max*y(5)-y(5)*(((phi_max-phi_min)*(conc(5)/MIC).^Hill_k)/((conc(5)/MIC).^Hill_k-(phi_min/(phi_max))));
        
        s(6)=phi_max*y(6)-y(6)*(((phi_max-phi_min)*(conc(6)/MIC).^Hill_k)/((conc(6)/MIC).^Hill_k-(phi_min/(phi_max))));
        
        s(7)=phi_max*y(7)-y(7)*(((phi_max-phi_min)*(conc(7)/MIC).^Hill_k)/((conc(7)/MIC).^Hill_k-(phi_min/(phi_max))));
        s(8)=phi_max*y(8)-y(8)*(((phi_max-phi_min)*(conc(8)/MIC).^Hill_k)/((conc(8)/MIC).^Hill_k-(phi_min/(phi_max))));
        s(9)=phi_max*y(9)-y(9)*(((phi_max-phi_min)*(conc(9)/MIC).^Hill_k)/((conc(9)/MIC).^Hill_k-(phi_min/(phi_max))));
       s(10)=phi_max*y(10)-y(10)*(((phi_max-phi_min)*(conc(10)/MIC).^Hill_k)/((conc(10)/MIC).^Hill_k-(phi_min/(phi_max))));
%         s(11)=phi_max*y(11)-y(11)*(((phi_max-phi_min)*(conc(11)/MIC).^Hill_k)/((conc(11)/MIC).^Hill_k-(phi_min/(phi_max))));
   
    end
 
output=values';

end
