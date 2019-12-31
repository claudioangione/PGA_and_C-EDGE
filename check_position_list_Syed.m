%% check position (in the multidimensional scaling) of the 180 cancer-related genes by Syed Haider
figure
load('cancer_genes_Syed.mat');
dist = zeros(numel(cancer_genes_Syed),1);
ix =[];
for i =1:numel(cancer_genes_Syed)
    new_ind = find(strcmp(geni_names,cancer_genes_Syed{i}));
    ix = [ix new_ind']; %list of positions in geni of the 96 cancer genes by Palsson (supplementary material file S1). Some genes may have multiple positions due to the transcriptional variants

end

C = distance_from_zero;  %colour according to the low edge score
scatter(Y(:,1),Y(:,2),200,C,'.');
title('Distance from (0,0)   [X = 96 cancer genes by Palsson]');
colorbar;
hold on
scatter(Y(ix,1),Y(ix,2),200,'black','X');