function [clearence_Times,gentaExtraVal,timeRunVal,total_BacteriaVal,combinationMIC]=dualTreatmentGentaAndAzithro(LHSSamples,clearedOutcomes_untreatedModel)

numLHSsamples=size(LHSSamples,1);
MICgentamicin=[4,8,16];
MICazithro=[0.5,1];

numMICgenta=length(MICgentamicin);
numMICazithro=length(MICazithro);
numComb=length(MICgentamicin)*length(MICazithro);%number of MIC combinations
combinationMIC=NaN(numComb,2);%at each iteration stores as a pair what the GEN and azitrho MICs are
clearence_Times=NaN(numLHSsamples,numComb);%colums are the MIC values and the rows gives the clearance times for the number of LHS samples
timeRunVal=cell(numComb,numLHSsamples); %rows are the MIC values and the columns are the LHS
total_BacteriaVal=cell(numComb,numLHSsamples); %rows are the MIC values and the columns are the LHS
gentaExtraVal=cell(numComb,numLHSsamples); 
initDose_1=240;%gentamicin dose.initial dose in mg. LHS samples generate Vd. The initial dose is fixed but the initial concnetration changes as the Vd changes.
initDose_2=1000;%azithromycin dose. LHS generates bioavailability and Vd.
iter=1;


for i=1:numMICgenta
    MICvalgenta=MICgentamicin(i);
    
    for k=1:numMICazithro
        MICvalazithro=MICazithro(k);
        combinationMIC(iter,1)=MICvalgenta;
        combinationMIC(iter,2)=MICvalazithro;
        for j=1:numLHSsamples
            paras=LHSSamples(j,:);
            paras(:,19)=MICvalgenta;%substituting gentamicin MIC
            paras(:,27)=MICvalazithro;%substituting for azithro MIC
            peakTime=clearedOutcomes_untreatedModel(j,3);
            MyPara=[paras,initDose_1,initDose_2];%adding MIC to the end
            [total_Bacteria,gentaExtra,clearenceTime_with_antibiotic,timeRun]=treatment_strategy_Genta_and_Azithro(MyPara,peakTime*24);
            clearence_Times(j,iter)=clearenceTime_with_antibiotic-peakTime;
            timeRunVal{iter,j}=timeRun;
            total_BacteriaVal{iter,j}=total_Bacteria;
            gentaExtraVal{iter,j}=gentaExtra;
        end
    iter=iter+1;    
    end

end
end
