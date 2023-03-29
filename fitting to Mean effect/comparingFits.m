%clear

%comparing the Hill function behaviour when fitted to only 1st experiment
%data and 2nd experiment data 

%phi_max-1st column.
%phi_min-2nd column
%MIC-3rd column
%Hill coefficient-4th column

para_Exp_1=[0.822952543	-9.999999761	0.193100246	1.668514688
0.661679061	-10	5.960922127	1.120426839
0.687082615	-2.05735811	0.028053751	0.910541624
0.658063209	-1.495249537	0.004088451	1.349808755
0.65996883	-0.53925113	0.00023388	1.746759206
0.780956815	-0.519424011	0.000123326	1.820136394
0.784260144	-0.40243278	0.352303814	0.893025235
0.732953718	-0.261938201	0.329852261	0.943958654
0.904036841	-9.999999324	0.001962493	0.933712553
];

para_Exp_2=[0.911860505	-6.34626587	0.302490069	1.762106731
0.848024201	-9.999999995	10	0.832291958
0.676967305	-0.98750439	0.024926933	0.973373293
0.917809933	-1.052251661	0.002758531	0.857338687
0.772949639	-0.364799809	0.00038365	1.76222666
0.670371237	-0.7618511	0.000462767	0.869094833
0.589742635	-0.068971268	0.998105864	1.641258996
0.799279323	-0.115615554	1.366059819	0.772413329
0.727484748	-11.54986306	0.001619076	0.900859407
];





conc=[0.000244141	0.004882813	3.05E-05	0.00012207	7.63E-06	1.95E-06	0.000488281	0.000244141	1.53E-05
0.002441406	0.048828125	0.000305176	0.001220703	7.63E-05	1.95E-05	0.004882813	0.002441406	0.000152588
0.004882813	0.09765625	0.000610352	0.002441406	0.000152588	3.91E-05	0.009765625	0.004882813	0.000305176
0.009765625	0.1953125	0.001220703	0.004882813	0.000305176	7.81E-05	0.01953125	0.009765625	0.000610352
0.01953125	0.390625	0.002441406	0.009765625	0.000610352	0.00015625	0.0390625	0.01953125	0.001220703
0.0390625	0.78125	0.004882813	0.01953125	0.001220703	0.0003125	0.078125	0.0390625	0.002441406
0.078125	1.5625	0.009765625	0.0390625	0.002441406	0.000625	0.15625	0.078125	0.004882813
0.15625	3.125	0.01953125	0.078125	0.004882813	0.00125	0.3125	0.15625	0.009765625
0.3125	6.25	0.0390625	0.15625	0.009765625	0.0025	0.625	0.3125	0.01953125
0.625	12.5	0.078125	0.3125	NaN	NaN	1.25	0.625	0.0390625
1.25	25	0.15625	NaN	NaN	NaN	2.5	1.25	0.078125
2.5	50	0.3125	NaN	NaN	NaN	NaN	2.5	0.15625
];%each coulmn gives the concnetrations of gentamicin	spectinomycin	azithromycin	penicillin	ceftriaxone	cefixime	chloro	tetra	cipro

numAntibiotics=9;
HillVals_exp_1=NaN(12,numAntibiotics);% coulmns are the antibioitcs. Rows give the Hill function value at each of the concentrations in conc
HillVals_exp_2=NaN(12,numAntibiotics);
HillVals_exp_both=NaN(12,numAntibiotics);

for i=1:numAntibiotics
    phi_max_1=para_Exp_1(i,1);
    phi_min_1=para_Exp_1(i,2);
    MIC_1=para_Exp_1(i,3);
    Hill_k_1=para_Exp_1(i,4);
    hill_function_exp_1=phi_max_1-((phi_max_1-phi_min_1).*(conc(:,i)/MIC_1).^Hill_k_1./((conc(:,i)/MIC_1).^Hill_k_1-(phi_min_1/phi_max_1)));
    numRows=length(hill_function_exp_1);
    if  numRows<13
       hill_function_exp_1( numRows+1:end)=NaN;
    end
    HillVals_exp_1(:,i)=hill_function_exp_1;
    
    
    phi_max_2=para_Exp_2(i,1);
    phi_min_2=para_Exp_2(i,2);
    MIC_2=para_Exp_2(i,3);
    Hill_k_2=para_Exp_2(i,4);
    hill_function_exp_2=phi_max_2-((phi_max_2-phi_min_2).*(conc(:,i)/MIC_2).^Hill_k_2./((conc(:,i)/MIC_2).^Hill_k_2-(phi_min_2/phi_max_2)));
    numRows=length(hill_function_exp_1);
    if  numRows<13
       hill_function_exp_2( numRows+1:end)=NaN;
    end
    HillVals_exp_2(:,i)=hill_function_exp_2;
       
      
end


