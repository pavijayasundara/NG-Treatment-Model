%load('cef1000azithro1000.mat')
samples_cleared_after_treatemnt=NaN(1,5);
for i=1:size(clearence_Times,2)
    indexes_cleared_after_treatemnt= find(clearence_Times(:,i)<7);
    samples_cleared_after_treatemnt(1,i)=(length(indexes_cleared_after_treatemnt)./5402)*100;
end
