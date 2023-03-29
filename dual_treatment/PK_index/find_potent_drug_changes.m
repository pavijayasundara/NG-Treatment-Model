function potent_indices = find_potent_drug_changes(potent_drug)

%%this gives the indexes in which the drug potency changes between GEN and
%%AZM

potentSample = potent_drug;
potent_indices = [];
tmp = 0;
while length(potentSample) >1
    
    if potentSample(1,1) == 1
        init_index = find(potentSample == 4, 1);
    end

    if potentSample(1,1) == 4
        init_index = find(potentSample == 1, 1);
    end
    potentSample = potentSample(init_index +1:end,1);%removing the indices before init_index
    potent_indices = [potent_indices; init_index+tmp];
    tmp = init_index+tmp;
end

end