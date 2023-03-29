function [unattached_NG,attached_NG,NGwithinPMN,NGwithinEpi,PMN,antibiotic,total_Bacteria,clearenceTime_with_antibiotic,timeRun]=treatment_strategy(LHSsamples,peakTime,initDose)

numLHSsamples=size(LHSsamples,1);
total_Bacteria=cell(1,numLHSsamples);
unattached_NG=cell(1,numLHSsamples);
attached_NG=cell(1,numLHSsamples);
NGwithinPMN=cell(1,numLHSsamples);
NGwithinEpi=cell(1,numLHSsamples);
PMN=cell(1,numLHSsamples);
antibiotic=cell(1,numLHSsamples);
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
    y0=[b,b_n,b_i,n,b_a,a];%all initial conditions
    tstart=0;
    
    timebetweenDose=24; %in hours
    second_doseTime=peakTime+timebetweenDose;
    third_doseTime=second_doseTime+timebetweenDose;
    fourth_doseTime=third_doseTime+timebetweenDose;
    fifth_doseTime=fourth_doseTime+timebetweenDose;
    sixth_doseTime=fifth_doseTime+timebetweenDose;
   
    %%% depending on the treatment strategy the number of doses that are
    %%% considered can be selected
    dosingTimes=[peakTime,5000];
    %dosingTimes=[peakTime,second_doseTime,2000];
    %dosingTimes=[peakTime,second_doseTime,third_doseTime,5000];
    %dosingTimes=[peakTime,second_doseTime,third_doseTime,fourth_doseTime,fifth_doseTime,sixth_doseTime,5000];
    options=odeset('Events',@eventfunction);
    tout=tstart;
    yout=y0;
    
    dose_1=initDose*paraValues(15)*0.76/(100*188.7);%first dose in micrograms per ml. Vd is fixed at 188.7. bd is generated from the LHS. free fraction 0.76
    dose_2=dose_1;
    dose_3=dose_1;  
    paras=paraValues;
    paras(:,15)=[];%removing the bd column and keeping all other parameters
    paras(:,18)=[];%removing the MIC column from LHS and keeping all other parameters
    step=0.015;%step size of Euler
    %Now the treatment para will be: alpha-15, phi_min-16,Hillk-17,delta-18, MIC-19
    for i=1:length(dosingTimes)
        tfinal=dosingTimes(i);
        [t,y,te]=eulerMethod(@GEP_treatment_model,tstart,y0,step,tfinal,paras);
        
        %%Accumulate output.
        nt=length(t);
        tout=[tout;t(2:nt)];
        yout=[yout;y(2:nt,:)];
        y0(1:5)=y(nt,1:5);
        if i==1
             y0(6)=y(nt,6)+dose_1;%adding dose 1 
        end

        if i==2
            y0(6)=y(nt,6)+dose_2;%adding 2nd dose
        end
        
        if i==3
            y0(6)=y(nt,6)+dose_3;%adding 3rd dose
        end
        if i==4
            y0(6)=y(nt,6)+dose_3;%adding 4th dose
        end
        if i==5
            y0(6)=y(nt,6)+dose_3;%adding 5th dose
        end
        if i==6
            y0(6)=y(nt,6)+dose_3;%adding 6th dose
        end
        
        tstart =t(nt);
        totNG=y(nt,1)+y(nt,2)+y(nt,3)+y(nt,5);
        if totNG<10
            break;
        end
        
    end
  
    total_NG=yout(:,1)+yout(:,2)+yout(:,3)+yout(:,5);
   
    total_Bacteria{1,j}=total_NG;
    if ~(isempty(te))
     clearenceTime_with_antibiotic{1,j}=te/24;
    end
    timeRun{1,j}=tout;
    unattached_NG{1,j}=yout(:,1);
    attached_NG{1,j}=yout(:,5);
    NGwithinPMN{1,j}=yout(:,2);
    NGwithinEpi{1,j}=yout(:,3);
    PMN{1,j}=yout(:,4);
    antibiotic{1,j}=yout(:,6);

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
    
    c=unattached_NG;
    maxlength=max(cellfun('size',c,1));
    for i=1:length(c)
        for j=cellfun('size',c(i),1)+1:maxlength
            c{i}(j)=NaN;
        end
    end
    unattached_NG=cell2mat(c);
    
    c=attached_NG;
    maxlength=max(cellfun('size',c,1));
    for i=1:length(c)
        for j=cellfun('size',c(i),1)+1:maxlength
            c{i}(j)=NaN;
        end
    end
    attached_NG=cell2mat(c);
    
    c=NGwithinPMN;
    maxlength=max(cellfun('size',c,1));
    for i=1:length(c)
        for j=cellfun('size',c(i),1)+1:maxlength
            c{i}(j)=NaN;
        end
    end
   NGwithinPMN=cell2mat(c);
 
   
    c=  NGwithinEpi;
    maxlength=max(cellfun('size',c,1));
    for i=1:length(c)
        for j=cellfun('size',c(i),1)+1:maxlength
            c{i}(j)=NaN;
        end
    end
    NGwithinEpi=cell2mat(c);
    
    c=PMN;
    maxlength=max(cellfun('size',c,1));
    for i=1:length(c)
        for j=cellfun('size',c(i),1)+1:maxlength
            c{i}(j)=NaN;
        end
    end
    PMN=cell2mat(c);
    
    c= antibiotic;
    maxlength=max(cellfun('size',c,1));
    for i=1:length(c)
        for j=cellfun('size',c(i),1)+1:maxlength
            c{i}(j)=NaN;
        end
    end
    antibiotic=cell2mat(c);
     
end
