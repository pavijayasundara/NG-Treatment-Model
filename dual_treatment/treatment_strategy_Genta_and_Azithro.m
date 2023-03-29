function [total_Bacteria,gentaExtra,clearenceTime_with_antibiotic,timeRun]=treatment_strategy_Genta_and_Azithro(LHSsamples,peakTime)

numLHSsamples=size(LHSsamples,1);
total_Bacteria=cell(1,numLHSsamples);
gentaExtra=cell(1,numLHSsamples);
clearenceTime_with_antibiotic=cell(1,numLHSsamples);
timeRun=cell(1,numLHSsamples);

for j=1:numLHSsamples
    
    b=1000;%initial size of the unattached bacteria
    b_n=0;%initial size of the bacteria within PMN
    b_i=0;%initial size of the bacteria within epithelial cells
    b_a=0;%initial size of the attached bacteria
    a=0;%antibiotic concentration (micro grams per mL)
    n=10^-8;%initial neutrophils
    paraValues=LHSsamples(j,:);
    y0=[b,b_n,b_i,n,b_a,a,0,0,0];%all initial conditions
    tstart=0;
    
    timebetweenDose=24;%in hours
    second_doseTime=peakTime+timebetweenDose;
    third_doseTime=second_doseTime+timebetweenDose;
    fourth_doseTime=third_doseTime+timebetweenDose;
    fifth_doseTime=fourth_doseTime+timebetweenDose;
    sixth_doseTime=fifth_doseTime+timebetweenDose;
    seventh_doseTime=sixth_doseTime+timebetweenDose;

    dosingTimes=[peakTime,5000];
    %dosingTimes=[peakTime,second_doseTime,5000];
    %dosingTimes=[peakTime,second_doseTime,third_doseTime,5000];
    %dosingTimes=[peakTime,second_doseTime,third_doseTime,fourth_doseTime,5000];
    %dosingTimes=[peakTime,second_doseTime,third_doseTime,fourth_doseTime,fifth_doseTime,sixth_doseTime,5000];
    %dosingTimes=[peakTime,second_doseTime,third_doseTime,fourth_doseTime,fifth_doseTime,sixth_doseTime,seventh_doseTime,5000];
    
    tout=tstart;
    yout=y0;
    
    dose_1=paraValues(33)/paraValues(15);%gentamicin dose in micrograms per ml
    dose_2=(paraValues(34)*0.37*0.88)/paraValues(16);%azithromycin dose in micrograms per ml
    
    paras=paraValues;
 
    parasAdjusted=paras([1:14 17:32]);%removing Vd genta and vd azithro
    step=0.05;%step size of Euler
    for i=1:length(dosingTimes)
        tfinal=dosingTimes(i);
        [t,y,te]=eulerMethod(@treatment_model_Loewe,tstart,y0,step,tfinal,parasAdjusted);
        %Accumulate output.
        nt=length(t);
        tout=[tout;t(2:nt)];
        yout=[yout;y(2:nt,:)];
        y0(1:5)=y(nt,1:5);
        y0(7)=y(nt,7);
        y0(8)=y(nt,8);
        if i==1
            y0(6)=y(nt,6)+dose_1;%adding gentamicin dose
            y0(9)=y(nt,9)+dose_2;%adding azithromycin dose
        end
        
        if i==2
            y0(6)=y(nt,6)+dose_1;%adding gentamicin dose
            y0(9)=y(nt,9)+dose_2;%azithromycin second dose
        end
        
        if i==3
            y0(6)=y(nt,6)+dose_1;%adding gentamicin dose
            y0(9)=y(nt,9);%azithromycin dose is not added
        end
        
        if i==4
            y0(6)=y(nt,6)+dose_1;%adding gentamicin dose
            y0(9)=y(nt,9);%azithromycin dose is not added
        end
        
        if i==5
            y0(6)=y(nt,6)+dose_1;%adding gentamicin dose
            y0(9)=y(nt,9);%azithromycin dose is not added
        end
        
        if i==6
            y0(6)=y(nt,6)+dose_1;%adding gentamicin dose
            y0(9)=y(nt,9);%azithromycin dose is not added
        end
        
        if i==7
            y0(6)=y(nt,6)+dose_1;%adding gentamicin dose
            y0(9)=y(nt,9);%azithromycin dose is not added
        end
        tstart =t(nt);
        
        totNG=y(nt,1)+y(nt,2)+y(nt,3)+y(nt,5);
        if totNG<10
            break;
        end
    end
    
    total_NG=yout(:,1)+yout(:,2)+yout(:,3)+yout(:,5);
    
    %numRows=length(total_NG);
    total_Bacteria{1,j}=total_NG;
    if ~(isempty(te))
        clearenceTime_with_antibiotic{1,j}=te/24;
    end
    timeRun{1,j}=tout;
    gentaExtra{1,j}=yout(:,6);
    
end


%converting cell to matrix
c=total_Bacteria;
maxlength=max(cellfun('size',c,1));
for i=1:length(c)
    for j=cellfun('size',c(i),1)+1:maxlength
        c{i}(j)=NaN;
    end
end
total_Bacteria=cell2mat(c);

c=gentaExtra;
maxlength=max(cellfun('size',c,1));
for i=1:length(c)
    for j=cellfun('size',c(i),1)+1:maxlength
        c{i}(j)=NaN;
    end
end
gentaExtra=cell2mat(c);

c=clearenceTime_with_antibiotic;
maxlength=max(cellfun('size',c,1));
for i=1:length(c)
    for j=cellfun('size',c(i),1)+1:maxlength
        c{i}(j)=NaN;
    end
end

clearenceTime_with_antibiotic=cell2mat(c);

c=timeRun;
maxlength=max(cellfun('size',c,1));
for i=1:length(c)
    for j=cellfun('size',c(i),1)+1:maxlength
        c{i}(j)=NaN;
    end
end
timeRun=cell2mat(c);
end


