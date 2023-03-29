function [gentaExtra,azithroExtra,gentaIntra, azithroIntra,timeRun, LoeweExtra, LoeweIntra, PotentExtra, PotentIntra]=treatment_strategy_Genta_and_Azithro(LHSsamples,peakTime)

numLHSsamples=size(LHSsamples,1);
gentaExtra=cell(1,numLHSsamples);
azithroExtra =cell(1,numLHSsamples);
gentaIntra=cell(1,numLHSsamples);
azithroIntra = cell(1,numLHSsamples);
timeRun=cell(1,numLHSsamples);
LoeweExtra = [];
LoeweIntra = [];
PotentExtra = [];
PotentIntra = [];

for j=1:numLHSsamples
    
    a=0;%antibiotic concentration (micro grams per mL)
    paraValues=LHSsamples(j,:);
    y0=[a,0,0,0];%all initial conditions
    tstart=0;
    
    timebetweenDose=8;%in hours
    second_doseTime=peakTime+timebetweenDose;
    third_doseTime=second_doseTime+timebetweenDose;
    fourth_doseTime=third_doseTime+timebetweenDose;
    fifth_doseTime=fourth_doseTime+timebetweenDose;
    sixth_doseTime=fifth_doseTime+timebetweenDose;
    seventh_doseTime=sixth_doseTime+timebetweenDose;
    
    dosingTimes=[peakTime,1000];
    %dosingTimes=[peakTime,second_doseTime,1000];
    %dosingTimes=[peakTime,second_doseTime,third_doseTime,1000];
    %dosingTimes=[peakTime,second_doseTime,third_doseTime,fourth_doseTime,5000];
    %dosingTimes=[peakTime,second_doseTime,third_doseTime,fourth_doseTime,fifth_doseTime,sixth_doseTime,5000];
    %dosingTimes=[peakTime,second_doseTime,third_doseTime,fourth_doseTime,fifth_doseTime,sixth_doseTime,seventh_doseTime,5000];
    
    
    tout=tstart;
    yout=y0;
    
    dose_1=paraValues(33)/paraValues(15);%gentamicin dose in micrograms per ml
    dose_2=(paraValues(34)*0.37*0.88)/paraValues(16);%azithromycin dose in micrograms per ml
    
    paras=paraValues;
    parasAdjusted=paras([1:14 17:32]);%removing Vd genta and vd azithro
    step=0.1;%step size of Euler
    
    for i=1:length(dosingTimes)
        tfinal=dosingTimes(i);
        [t,y,~, ~, ~]=eulerMethod(@treatment_model_Loewe,tstart,y0,step,tfinal,parasAdjusted);
        %Accumulate output.
        nt=length(t);
        tout=[tout;t(2:nt)];
        yout=[yout;y(2:nt,:)];
        
        Loewe_conc_extra = zeros(nt, 1);
        Loewe_conc_intra = zeros(nt, 1);
        potent_ex = zeros(nt, 1);
        potent_in = zeros(nt, 1);
        
        
        for K = 1 : nt
            [~, Loewe_conc_extra(K), Loewe_conc_intra(K), potent_ex(K), potent_in(K)] = treatment_model_Loewe(t(K), y(K,:),parasAdjusted);
        end
        LoeweExtra = [LoeweExtra; Loewe_conc_extra];
        LoeweIntra = [LoeweIntra; Loewe_conc_intra];
        PotentExtra = [PotentExtra; potent_ex];
        PotentIntra = [PotentIntra; potent_ex];
        
        y0(2)=y(nt,2);
        y0(3)=y(nt,3);
        
        if i==1
            y0(1)=y(nt,1)+dose_1;%adding gentamicin dose
            y0(4)=y(nt,4)+dose_2;%adding azithromycin dose
        end
        
        if i==2
            y0(1)=y(nt,1)+dose_1;%adding gentamicin dose
            y0(4)=y(nt,4)+dose_2;%azithromycin second dose
        end
        
        if i==3
            y0(1)=y(nt,1)+dose_1;%adding gentamicin dose
            y0(4)=y(nt,4);%azithromycin dose is not added
        end
        
        if i==4
            y0(1)=y(nt,1)+dose_1;%adding gentamicin dose
            y0(4)=y(nt,4);%azithromycin dose is not added
        end
        
        if i==5
            y0(1)=y(nt,1)+dose_1;%adding gentamicin dose
            y0(4)=y(nt,4);%azithromycin dose is not added
        end
        
        if i==6
            y0(1)=y(nt,1)+dose_1;%adding gentamicin dose
            y0(4)=y(nt,4);%azithromycin dose is not added
        end
        
        if i==7
            y0(1)=y(nt,1)+dose_1;%adding gentamicin dose
            y0(4)=y(nt,4);%azithromycin dose is not added
        end
        tstart =t(nt);
        
    end
    
    timeRun{1,j}=tout;
    gentaExtra{1,j}=yout(:,1);
    gentaIntra{1,j}=yout(:,2);
    azithroExtra{1,j}=yout(:,4);
    azithroIntra{1,j}=yout(:,3);
end
LoeweExtra(end)=[];
LoeweIntra(end)= [];

PotentExtra(end)=[];
PotentIntra(end)= [];


c=gentaExtra;
maxlength=max(cellfun('size',c,1));
for i=1:length(c)
    for j=cellfun('size',c(i),1)+1:maxlength
        c{i}(j)=NaN;
    end
end
gentaExtra=cell2mat(c);

c=gentaIntra;
maxlength=max(cellfun('size',c,1));
for i=1:length(c)
    for j=cellfun('size',c(i),1)+1:maxlength
        c{i}(j)=NaN;
    end
end
gentaIntra=cell2mat(c);

c=azithroExtra;
maxlength=max(cellfun('size',c,1));
for i=1:length(c)
    for j=cellfun('size',c(i),1)+1:maxlength
        c{i}(j)=NaN;
    end
end
azithroExtra=cell2mat(c);

c=azithroIntra;
maxlength=max(cellfun('size',c,1));
for i=1:length(c)
    for j=cellfun('size',c(i),1)+1:maxlength
        c{i}(j)=NaN;
    end
end
azithroIntra=cell2mat(c);

c=timeRun;
maxlength=max(cellfun('size',c,1));
for i=1:length(c)
    for j=cellfun('size',c(i),1)+1:maxlength
        c{i}(j)=NaN;
    end
end
timeRun=cell2mat(c);

end


