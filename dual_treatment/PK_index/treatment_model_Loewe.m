    function [s, c1, c2, ex_potent, in_potent]= treatment_model_Loewe(~,y,para)

    %Dual treatemnt with gentamicin and azithromycin. 
    r1=para(1);%replication rate  of unattached bacteria

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
    

    s=zeros(4,1);
    
    s(1)=-k_10_1*y(1)-k_12_1*y(1)+k_21_1*y(2)*(V_c1/V_e1);%gentamicin extracellular
    s(2)=k_12_1*y(1)*(V_e1/V_c1)-k_21_1*y(2);%gentamicin intracellular

    s(3)=k_12_2*y(4)*(V_e2/V_c2)-k_21_2*y(3);%azithromycin intracellular
    s(4)=-k_10_2*y(4)-k_12_2*y(4)+k_21_2*y(3)*(V_c2/V_e2); %azithromycin extracellular

    [hill_function_gentaEx,hill_function_azithroEx]=hillEffect(y(1),y(4),parasGenta,parasAzithro);
    [hill_function_gentaIn,hill_function_azithroIn]=hillEffect(y(2),y(3),parasGenta,parasAzithro);

    %finding the highly potent drug in the extracellular space.
    if hill_function_gentaEx>hill_function_azithroEx%here gentamicin is the more potent drug extracellularly

        effect_azithro_extra=(r1-phi_min_2)*(y(4)/MIC_2)^Hillk_2/((y(4)/MIC_2)^Hillk_2-(phi_min_2/(r1)));
        c_eq= MIC_1*(phi_min_1*effect_azithro_extra/(r1*(effect_azithro_extra-(r1-phi_min_1))))^(1/Hillk_1);
        c1=y(1)+c_eq;
        ex_potent = 4; % because genta is more potent the MIC we should look at is Genta which is 4

    else %here azithromycin is the highly potent drug extracellularly
        effect_genta_extra=(r1-phi_min_1)*(y(1)/MIC_1)^Hillk_1/((y(1)/MIC_1)^Hillk_1-(phi_min_1/(r1)));
        c_eq= MIC_2*(phi_min_2*effect_genta_extra/(r1*(effect_genta_extra-(r1-phi_min_2))))^(1/Hillk_2);
        c1=y(4)+c_eq;
        ex_potent = 1; % because AZM is more potent the MIC we should look at is AZM which is 1


    end

    %finding the highly potent drug in the intracellular space.
    if hill_function_gentaIn>hill_function_azithroIn%here gentamicin is the more potent drug intracellularly

        effect_azithro_intra=(r1-phi_min_2)*(y(3)/MIC_2)^Hillk_2/((y(3)/MIC_2)^Hillk_2-(phi_min_2/(r1)));
        c_eq= MIC_1*(phi_min_1*effect_azithro_intra/(r1*(effect_azithro_intra-(r1-phi_min_1))))^(1/Hillk_1);
        c2=y(2)+c_eq;
        in_potent = 4; % because genta is more potent the MIC we should look at is Genta which is 4

    else %here azithromycin is the highly potent drug intracellularly
        effect_genta_intra=(r1-phi_min_1)*(y(2)/MIC_1)^Hillk_1/((y(2)/MIC_1)^Hillk_1-(phi_min_1/(r1)));
        c_eq= MIC_2*(phi_min_2*effect_genta_intra/(r1*(effect_genta_intra-(r1-phi_min_2))))^(1/Hillk_2);
        c2=y(3)+c_eq;
        in_potent = 1; % because AZM is more potent the MIC we should look at is AZM which is 1
    end

    end