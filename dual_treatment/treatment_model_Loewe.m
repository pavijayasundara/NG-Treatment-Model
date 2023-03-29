    function s= treatment_model_Loewe(~,y,para)
    %Daul treatemnt with gentamicin and azithromycin. 

    r1=para(1);%replication rate  of unattached bacteria
    a1=para(4);%bacterial attachment rate to epithelial
    p=para(11);%proportion surving within PMN
    d=para(9);%rate that NG are engulfed by PMN
    d2=para(12);%washout rate of unattached bacteria
    N_max=2.5*10^10;%maximum number of neutrophils in the body
    mu=para(7);%neutrophil activation rate
    d3=para(13);%death rate of neutrophil
    k=10^7;%carrying capacity
    a2=1/para(5);%infection rate of epithelial cells
    k2=k/a2;%carrying capacity of epithelial cell
    c=para(10);%ratio dependent constant
    e=para(6);%rate of exit from epithelial cells
    k3=para(14);%carrying capacity of PMN
    r2=para(2);%replication rate of internalized NG within epithelial cells
    r3=para(3);%replication of survinvg NG within PMN
    eta=para(8);%proportion of NG internalized

    %gentamicin parameters
    phi_min_1=para(15);
    Hillk_1=para(16);%Hill coefficient
    MIC_1=para(17);
    k_10_1=para(18);
    k_12_1=para(19);%transfer form extracellualr to intracellular rate constant
    k_21_1=para(20);%transfer form intracellualr to extracellular rate constant
    V_e1=para(21);
    V_c1=para(22);

    %azitromycin parameters
    phi_min_2=para(23);
    Hillk_2=para(24);%Hill coefficient
    MIC_2=para(25);
    k_10_2=para(26);%elimination rate
    k_12_2=para(27);%transfer form extracellualr to intracellular rate constant
    k_21_2=para(28);%transfer form intracellualr to extracellular rate constant
    V_e2=para(29);
    V_c2=para(30);
    parasGenta=[phi_min_1,Hillk_1,MIC_1];
    parasAzithro=[phi_min_2,Hillk_2,MIC_2];
    

    s=zeros(9,1);

    s(6)=-k_10_1*y(6)-k_12_1*y(6)+k_21_1*y(7)*(V_c1/V_e1);%gentamicin extracellular
    s(7)=k_12_1*y(6)*(V_e1/V_c1)-k_21_1*y(7);%gentamicin intracellular

    s(8)=k_12_2*y(9)*(V_e2/V_c2)-k_21_2*y(8);%azithromycn intracellular
    s(9)=-k_10_2*y(9)-k_12_2*y(9)+k_21_2*y(8)*(V_c2/V_e2); %azithromycn extracellular

    [hill_funciton_gentaEx,hill_funciton_azithroEx]=hillEffect(y(6),y(9),parasGenta,parasAzithro);
    [hill_funciton_gentaIn,hill_funciton_azithroIn]=hillEffect(y(7),y(8),parasGenta,parasAzithro);

    %finding the highly potent drug in the extracellular space.
    if hill_funciton_gentaEx>hill_funciton_azithroEx%here gentamicin is the more potent drug extracellularly

        effect_azithro_extra=(r1-phi_min_2)*(y(9)/MIC_2)^Hillk_2/((y(9)/MIC_2)^Hillk_2-(phi_min_2/(r1)));
        c_eq= MIC_1*(phi_min_1*effect_azithro_extra/(r1*(effect_azithro_extra-(r1-phi_min_1))))^(1/Hillk_1);
        c1=y(6)+c_eq;
        Emax_extra=(((r1-phi_min_1)*(c1/MIC_1)^Hillk_1)/((c1/MIC_1)^Hillk_1-(phi_min_1/(r1))));
    else %here azithromycin is the highly potent drug extracellularly
        effect_genta_extra=(r1-phi_min_1)*(y(6)/MIC_1)^Hillk_1/((y(6)/MIC_1)^Hillk_1-(phi_min_1/(r1)));
        c_eq= MIC_2*(phi_min_2*effect_genta_extra/(r1*(effect_genta_extra-(r1-phi_min_2))))^(1/Hillk_2);
        c1=y(9)+c_eq;
        
        Emax_extra=(((r1-phi_min_2)*(c1/MIC_2)^Hillk_2)/((c1/MIC_2)^Hillk_2-(phi_min_2/(r1))));

    end

    %finding the highly potent drug in the intracellular space.
    if hill_funciton_gentaIn>hill_funciton_azithroIn%here gentamicin is the more potent drug intracellularly

        effect_azithro_intra=(r1-phi_min_2)*(y(8)/MIC_2)^Hillk_2/((y(8)/MIC_2)^Hillk_2-(phi_min_2/(r1)));
        c_eq= MIC_1*(phi_min_1*effect_azithro_intra/(r1*(effect_azithro_intra-(r1-phi_min_1))))^(1/Hillk_1);
        c2=y(7)+c_eq;
        Emax_intra_PMN=(((r3-phi_min_1)*((c2)/MIC_1)^Hillk_1)/(((c2)/MIC_1)^Hillk_1-(phi_min_1/r3)));
        Emax_intra_epi=(((r2-phi_min_1)*((c2)/MIC_1)^Hillk_1)/(((c2)/MIC_1)^Hillk_1-(phi_min_1/r2)));
    else %here azithromycin is the highly potent drug intracellularly
        effect_genta_intra=(r1-phi_min_1)*(y(7)/MIC_1)^Hillk_1/((y(7)/MIC_1)^Hillk_1-(phi_min_1/(r1)));
        c_eq= MIC_2*(phi_min_2*effect_genta_intra/(r1*(effect_genta_intra-(r1-phi_min_2))))^(1/Hillk_2);
        c2=y(8)+c_eq;
     
        Emax_intra_PMN=(((r3-phi_min_2)*((c2)/MIC_2)^Hillk_2)/(((c2)/MIC_2)^Hillk_2-(phi_min_2/r3)));
        Emax_intra_epi=(((r2-phi_min_2)*((c2)/MIC_2)^Hillk_2)/(((c2)/MIC_2)^Hillk_2-(phi_min_2/r2)));
    end

    s(1)=(1-((y(1)+y(5))/k))*(r1*y(1)+d3*y(2)+e*y(3))-(d*y(1)*y(4)/(c*y(4))+y(1))-d2*y(1)...
        -a1*y(1)*(1-(y(5)/(k2)))....
        - y(1)*(1-((y(1)+y(5))/k))*Emax_extra;%rate of change of unattached bacteria
    s(2)=(1-(y(2)/(y(4)*k3)))*((p*d*y(4)*y(1)/(c*y(4)+y(1)))+(p*d*y(4)*y(5)/(c*y(4)+y(5)))+...
        r3*y(2))-d3*y(2)-....
        y(2)*(1-(y(2)/(y(4)*k3)))*Emax_intra_PMN;%NG within PMN
    s(3)=(1-(y(3)/k2))*(eta*y(5)+r2*y(3))-e*y(3)-....
        y(3)*(1-(y(3)/k2))*Emax_intra_epi;%NG within epithelial
    s(4)=mu*(N_max-y(4))*(y(1)+y(5))-d3*y(4);%activated neutrophil
    s(5)=(r1*y(5))*(1-((y(1)+y(5))/k))-(d*y(4)*y(5)/(c*y(4)+y(5)))-...
        eta*y(5)+a1*y(1)*(1-(y(5)/k2))....
        - y(5)*(1-((y(1)+y(5))/k))*Emax_extra;%attached bacteria
    end