function [gentaExtraVal,gentaIntraVal,azithroExtraVal, azithroIntraVal,timeRunVal,combinationMIC, Loewe_ExtraVal, Loewe_IntraVal, potent_ex_val, potent_in_val]=dualTreatmentGentaAndAzithro(LHSSamples,clearedOutcomes_untreatedModel)


numLHSsamples=size(LHSSamples,1);
MICgentamicin= 4;
MICazithro= 1;

numMICgenta=length(MICgentamicin);
numMICazithro=length(MICazithro);
numComb=length(MICgentamicin)*length(MICazithro);%number of MIC combinations
combinationMIC=NaN(numComb,2);%at each iteration stores as a pair what the GEN and AZM MICs are

timeRunVal=cell(numComb,numLHSsamples); 
gentaExtraVal=cell(numComb,numLHSsamples);
azithroExtraVal=cell(numComb,numLHSsamples);
gentaIntraVal=cell(numComb,numLHSsamples);
azithroIntraVal =cell(numComb,numLHSsamples);

Loewe_ExtraVal= cell(numComb,numLHSsamples); 
Loewe_IntraVal= cell(numComb,numLHSsamples); 

potent_ex_val= cell(numComb,numLHSsamples); 
potent_in_val= cell(numComb,numLHSsamples); 

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
            [gentaExtra,azithroExtra,gentaIntra, azithroIntra,timeRun, LoeweExtra, LoeweIntra,...
                potent_ex, potent_in]=treatment_strategy_Genta_and_Azithro(MyPara,peakTime*24);
           
            timeRunVal{iter,j}=timeRun;
            gentaExtraVal{iter,j}=gentaExtra;
            gentaIntraVal{iter,j}=gentaIntra;
            azithroExtraVal{iter,j}=azithroExtra;
            azithroIntraVal{iter,j}=azithroIntra;
            Loewe_ExtraVal{iter,j}=LoeweExtra;
            Loewe_IntraVal{iter,j}=LoeweIntra; 
            potent_ex_val{iter,j}=potent_ex;
            potent_in_val{iter,j}=potent_in; 
        end
    iter=iter+1;    
    end

end

timeRunVal = cell2mat(timeRunVal);
gentaExtraVal =cell2mat(gentaExtraVal);
gentaIntraVal =cell2mat(gentaIntraVal);
azithroExtraVal =cell2mat(azithroExtraVal);
azithroIntraVal =cell2mat(azithroIntraVal);
Loewe_ExtraVal =cell2mat(Loewe_ExtraVal);
Loewe_IntraVal =cell2mat(Loewe_IntraVal);
potent_ex_val = cell2mat(potent_ex_val);
potent_in_val = cell2mat(potent_in_val);
end
