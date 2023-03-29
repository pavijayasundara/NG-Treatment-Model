%the effects (from the Hill function) were estimated for individual fits to
%experiments 1 and 2.
%Then the mean "effect" was calculated from the two experiments.
%to these mean effect that varies over concentrations the Hill function was
%fitted by taking the effect-concentration values as the data points.
close all

clear
clc

comparingFits%script which evaluates the individual experiment Hill function

%Hill function effects from experiment 1
effect_experiment_1=HillVals_exp_1;
%Hill function effects from experiment 2

effect_experiment_2=HillVals_exp_2;
mean_effect=(effect_experiment_1+effect_experiment_2)/2;


titles={'gentamicin','spectinomycin','azithromycin','penicillin','ceftriaxone','cefixime','chloramphenicol','tetracycline','ciprofloxacin'};
numAntibioitcs=9;
MICest=NaN(numAntibioitcs,1);
Hill_k=NaN(numAntibioitcs,1);
errors=NaN(numAntibioitcs,1);
pi_min=NaN(numAntibioitcs,1);
pi_max=NaN(numAntibioitcs,1);


figure('units','normalized','outerposition',[0 0 1 1])

for i=1:9
   subplot(3,3,i)
    h1=plot(conc(:,i), effect_experiment_1(:,i),'b','linewidth',2);
    hold on
    h2=plot(conc(:,i), effect_experiment_2(:,i),'color',[255,140,0]./255,'linewidth',2,'linestyle','--');
 
   
    xlabel('conc. (mg/L)')
    %ylabel('net growth rate')
    set(gca,'fontsize',16,'fontname','Times','Box','off')
    if i==4
        [fittedVals,errorVals,confInt,h3]=fitting_mean_effect2(mean_effect(1:10,i),conc(1:10,i),titles{i},-8);
    elseif i==5
        [fittedVals,errorVals,confInt,h3]=fitting_mean_effect2(mean_effect(1:9,i),conc(1:9,i),titles{i},-8);   
    elseif i==6
        [fittedVals,errorVals,confInt,h3]=fitting_mean_effect2(mean_effect(1:9,i),conc(1:9,i),titles{i},-8);   
    elseif i==7
        [fittedVals,errorVals,confInt,h3]=fitting_mean_effect2(mean_effect(1:11,i),conc(1:11,i),titles{i},-8);   
    else
        [fittedVals,errorVals,confInt,h3]=fitting_mean_effect2(mean_effect(:,i),conc(:,i),titles{i},-8);
    end
    
     if i==1
       legend([h1,h2,h3],{'exp 1','exp 2','mean fit'},'Box','off') 
    end
%parameters estimated: phi-max,MIC and Hill_k,phi_min
    pi_max(i)=fittedVals(1);
    MICest(i)=fittedVals(2);
    Hill_k(i)=fittedVals(3);
    pi_min(i)=fittedVals(4);

    
end
results=[pi_max,pi_min,MICest,Hill_k];

annotation('textbox',...
    [0.090311859443631 0.90937019969278 0.0268184480234261 0.0430107526881719],...
    'String',{'(a)'},...
    'LineStyle','none',...
    'FontSize',16,...
    'FontName','Times',...
    'FitBoxToText','off');

% Create textbox
annotation('textbox',...
    [0.372156661786237 0.920122887864823 0.0268184480234261 0.0430107526881719],...
    'String','(b)',...
    'LineStyle','none',...
    'FontSize',16,...
    'FontName','Times',...
    'FitBoxToText','off');

% Create textbox
annotation('textbox',...
    [0.648144948755491 0.926267281105991 0.026818448023426 0.0430107526881719],...
    'String','(c)',...
    'LineStyle','none',...
    'FontSize',16,...
    'FontName','Times',...
    'FitBoxToText','off');

% Create textbox
annotation('textbox',...
    [0.0932401171303075 0.629800307219662 0.0268184480234261 0.0430107526881719],...
    'String','(d)',...
    'LineStyle','none',...
    'FontSize',16,...
    'FontName','Times',...
    'FitBoxToText','off');

% Create textbox
annotation('textbox',...
    [0.374352855051245 0.623655913978494 0.0268184480234263 0.0430107526881719],...
    'String','(e)',...
    'LineStyle','none',...
    'FontSize',16,...
    'FontName','Times',...
    'FitBoxToText','off');

% Create textbox
annotation('textbox',...
    [0.655465592972182 0.622119815668203 0.026818448023426 0.0430107526881719],...
    'String','(f)',...
    'LineStyle','none',...
    'FontSize',16,...
    'FontName','Times',...
    'FitBoxToText','off');

% Create textbox
annotation('textbox',...
    [0.0932401171303075 0.33026113671275 0.0268184480234261 0.0430107526881718],...
    'String','(g)',...
    'LineStyle','none',...
    'FontSize',16,...
    'FontName','Times',...
    'FitBoxToText','off');

% Create textbox
annotation('textbox',...
    [0.373620790629576 0.331797235023041 0.0268184480234263 0.0430107526881718],...
    'String','(h)',...
    'LineStyle','none',...
    'FontSize',16,...
    'FontName','Times',...
    'FitBoxToText','off');

% Create textbox
annotation('textbox',...
    [0.656197657393851 0.319508448540707 0.026818448023426 0.0430107526881719],...
    'String','(i)',...
    'LineStyle','none',...
    'FontSize',16,...
    'FontName','Times',...
    'FitBoxToText','off');

annotation('textarrow',[0.00541786930256006 0.00319564708033785],...
    [0.689516326837412 0.477664474985561],'TextRotation',90,...
    'String',{'net growth rate (h^{-1})'},...
    'LineStyle','none',...
    'HeadStyle','none',...
    'FontSize',18,...
    'FontName','Times');
set(gca,'fontsize',16,'fontname','Times','Box','off')