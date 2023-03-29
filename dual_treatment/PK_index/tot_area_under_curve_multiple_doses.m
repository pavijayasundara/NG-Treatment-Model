function totareaAboveMIC=tot_area_under_curve_multiple_doses(antibiotic, potent_drug)

numSamples=size(antibiotic,2);
totareaAboveMIC=NaN(1,numSamples);

for j=1:numSamples
    potent_indices = find_potent_drug_changes(potent_drug(:,j));
    potent_indices=[1;potent_indices];
    
    cumulativeArea=0;
    
    
    if length(potent_indices)== 1
        totalAreaBelow =trapz(antibiotic(:,j));
        potent_MIC = potent_drug(potent_indices(1),j);
        totareaAboveMIC(1,j)= totalAreaBelow/potent_MIC;
    else
        for k =1:length(potent_indices)-1

            totalAreaBelow =trapz(antibiotic(potent_indices(k):potent_indices(k+1)-1,j));
            potent_MIC = potent_drug(potent_indices(k),j);
            area_below_MIC = totalAreaBelow/potent_MIC ;
            cumulativeArea= cumulativeArea + area_below_MIC ;
        end
        totareaAboveMIC(1,j)= cumulativeArea;
        
    end

end

end